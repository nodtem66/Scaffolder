// Copyright (C) 2017, Dominik Wodniok
// This software may be modified and distributed under the terms
// of the BSD 3-Clause license.
// See the LICENSE.txt file for details.

#ifndef DUALMC_H_INCLUDED
#define DUALMC_H_INCLUDED

/// \file   dualmc.h
/// \author Dominik Wodniok
/// \date   2009

// c includes
#include <cstdint>

// stl includes
#include <unordered_map>
#include <vector>
#include "ProgressBar.hpp"
#include "oneapi/tbb.h"
#ifndef PROGRESS_BAR_COLUMN
#define PROGRESS_BAR_COLUMN 40
#endif
#define USE_PARALLEL

namespace dualmc {


	typedef double VertexComponentsType;
	typedef int32_t QuadIndexType;

	/// vertex structure for dual points
	struct Vertex {
		/// non-initializing constructor
		Vertex();

		/// initializing constructor
		Vertex(VertexComponentsType x, VertexComponentsType y, VertexComponentsType z);

		/// initializing constructor
		Vertex(Vertex const& v);

		// components
		VertexComponentsType x, y, z;
	};

	/// quad indices structure
	struct Quad {
		/// non-initializing constructor
		Quad();

		/// initializing constructor
		Quad(QuadIndexType i0, QuadIndexType i1, QuadIndexType i2, QuadIndexType i3);

		// quad indices
		QuadIndexType i0, i1, i2, i3;

	};

	/// \class  DualMC
	/// \author Dominik Wodniok
	/// \date   2009
	/// Class which implements the dual marching cubes algorithm from Gregory M. Nielson.
	/// Faces and vertices of the standard marching cubes algorithm correspond to
	/// vertices and faces in the dual algorithm. As a vertex in standard marching cubes
	/// usually is shared by 4 faces, the dual mesh is entirely made from quadrangles.
	/// Unfortunately, under rare circumstances the original algorithm can create
	/// non-manifold meshes. See the remarks of the original paper on this.
	/// The class optionally can guarantee manifold meshes by taking the Manifold
	/// Dual Marching Cubes approach from Rephael Wenger as described in
	/// chapter 3.3.5 of his book "Isosurfaces: Geometry, Topology, and Algorithms".
	template<class T> class DualMC {

	public:
		// typedefs
		typedef T VolumeDataType;

		// Callback function to send the progression
		std::function<void(int)> callback = NULL;

		/// Extracts the iso surface for a given volume and iso value.
		/// Output is a list of vertices and a list of indices, which connect
		/// vertices to quads.
		/// The quad mesh either uses shared vertex indices or is a quad soup if
		/// desired.
		void build(
			VolumeDataType const* data,
			int32_t const dimX, int32_t const dimY, int32_t const dimZ,
			VolumeDataType const iso,
			bool const generateManifold,
			bool const generateSoup,
			std::vector<Vertex>& vertices,
			std::vector<Quad>& quads
		);

	private:

		/// Extract quad mesh with shared vertex indices.
		void buildSharedVerticesQuads(
			VolumeDataType const iso,
			std::vector<Vertex>& vertices,
			std::vector<Quad>& quads
		);

		/// Extract quad soup.
		void buildQuadSoup(
			VolumeDataType const iso,
			std::vector<Vertex>& vertices,
			std::vector<Quad>& quads
		);


	private:

		/// enum with edge codes for a 12-bit voxel edge mask to indicate
		/// grid edges which intersect the ISO surface of classic marching cubes
		enum DMCEdgeCode {
			EDGE0 = 1,
			EDGE1 = 1 << 1,
			EDGE2 = 1 << 2,
			EDGE3 = 1 << 3,
			EDGE4 = 1 << 4,
			EDGE5 = 1 << 5,
			EDGE6 = 1 << 6,
			EDGE7 = 1 << 7,
			EDGE8 = 1 << 8,
			EDGE9 = 1 << 9,
			EDGE10 = 1 << 10,
			EDGE11 = 1 << 11,
			FORCE_32BIT = 0xffffffff
		};

		/// get the 8-bit in-out mask for the voxel corners of the cell cube at (cx,cy,cz)
		/// and the given iso value
		int getCellCode(int32_t const cx, int32_t const cy, int32_t const cz, VolumeDataType const iso) const;

		/// Get the 12-bit dual point code mask, which encodes the traditional
		/// marching cube vertices of the traditional marching cubes face which
		/// corresponds to the dual point.
		/// This is also where the manifold dual marching cubes algorithm is
		/// implemented.
		int getDualPointCode(int32_t const cx, int32_t const cy, int32_t const cz,
			VolumeDataType const iso, DMCEdgeCode const edge) const;

		/// Given a dual point code and iso value, compute the dual point.
		void calculateDualPoint(int32_t const cx, int32_t const cy, int32_t const cz,
			VolumeDataType const iso, int const pointCode, Vertex& v) const;

		/// Compute a linearized cell cube index.
		int32_t gA(int32_t const x, int32_t const y, int32_t const z) const;

	private:
		// static lookup tables needed for (manifold) dual marching cubes

		/// Dual Marching Cubes table
		/// Encodes the edge vertices for the 256 marching cubes cases.
		/// A marching cube case produces up to four faces and ,thus, up to four
		/// dual points.
		static int32_t const dualPointsList[256][4];

		/// Table which encodes the ambiguous face of cube configurations, which
		/// can cause non-manifold meshes.
		/// Needed for manifold dual marching cubes.
		static uint8_t const problematicConfigs[256];

	private:

		/// convenience volume extent array for x-,y-, and z-dimension
		int32_t dims[3];

		/// convenience volume data point
		VolumeDataType const* data;

		/// store whether the manifold dual marching cubes algorithm should be
		/// applied.
		bool generateManifold;

		/// Dual point key structure for hashing of shared vertices
		struct DualPointKey {
			// a dual point can be uniquely identified by ite linearized volume cell
			// id and point code
			int32_t linearizedCellID;
			int pointCode;
			/// Equal operator for unordered map
			bool operator==(DualPointKey const& other) const;
		};

		/// Functor for dual point key hash generation
		struct DualPointKeyHash {
			size_t operator()(DualPointKey const& k) const {
				return size_t(k.linearizedCellID) | (size_t(k.pointCode) << 32u);
			}
		};

		// TBB Mutex type
		typedef tbb::queuing_mutex MutexType;
		typedef tbb::queuing_mutex CallbackMutexType;


		/// Hash map for shared vertex index computations
		typedef std::unordered_map<DualPointKey, QuadIndexType, DualPointKeyHash> PointToIndexMap;
		//std::unordered_map<DualPointKey, QuadIndexType, DualPointKeyHash> pointToIndex;

		/// Get the shared index of a dual point which is uniquly identified by its
		/// cell cube index and a cube edge. The dual point is computed,
		/// if it has not been computed before.
		QuadIndexType getSharedDualPointIndex(int32_t const cx, int32_t const cy, int32_t const cz,
			VolumeDataType const iso, DMCEdgeCode const edge,
			std::vector<Vertex>& vertices, PointToIndexMap& p);
	};

	// inline function definitions

	//------------------------------------------------------------------------------

	inline
		Vertex::Vertex() {}

	//------------------------------------------------------------------------------

	inline
		Vertex::Vertex(
			VertexComponentsType x,
			VertexComponentsType y,
			VertexComponentsType z
		) : x(x), y(y), z(z) {}

	//------------------------------------------------------------------------------

	inline
		Vertex::Vertex(Vertex const& v) : x(v.x), y(v.y), z(v.z) {}

	//------------------------------------------------------------------------------

	inline
		Quad::Quad() {}

	//------------------------------------------------------------------------------

	inline
		Quad::Quad(
			QuadIndexType i0,
			QuadIndexType i1,
			QuadIndexType i2,
			QuadIndexType i3
		) : i0(i0), i1(i1), i2(i2), i3(i3) {}

	//------------------------------------------------------------------------------

	template<class T> inline
		int32_t DualMC<T>::gA(int32_t const x, int32_t const y, int32_t const z) const {
		return x + dims[0] * (y + dims[1] * z);
	}

	//------------------------------------------------------------------------------
	template<class T> inline
		bool DualMC<T>::DualPointKey::operator==(typename DualMC<T>::DualPointKey const& other) const {
		return linearizedCellID == other.linearizedCellID && pointCode == other.pointCode;
	}

	template<class T> inline
		int DualMC<T>::getCellCode(int32_t const cx, int32_t const cy, int32_t const cz, VolumeDataType const iso) const {
		// determine for each cube corner if it is outside or inside
		int code = 0;
		if (data[gA(cx, cy, cz)] >= iso)
			code |= 1;
		if (data[gA(cx + 1, cy, cz)] >= iso)
			code |= 2;
		if (data[gA(cx, cy + 1, cz)] >= iso)
			code |= 4;
		if (data[gA(cx + 1, cy + 1, cz)] >= iso)
			code |= 8;
		if (data[gA(cx, cy, cz + 1)] >= iso)
			code |= 16;
		if (data[gA(cx + 1, cy, cz + 1)] >= iso)
			code |= 32;
		if (data[gA(cx, cy + 1, cz + 1)] >= iso)
			code |= 64;
		if (data[gA(cx + 1, cy + 1, cz + 1)] >= iso)
			code |= 128;
		return code;
	}

	//------------------------------------------------------------------------------

	template<class T> inline
		int DualMC<T>::getDualPointCode(int32_t const cx, int32_t const cy, int32_t const cz, VolumeDataType const iso, DMCEdgeCode const edge) const {
		int cubeCode = getCellCode(cx, cy, cz, iso);

		// is manifold dual marching cubes desired?
		if (generateManifold) {
			// The Manifold Dual Marching Cubes approach from Rephael Wenger as described in
			// chapter 3.3.5 of his book "Isosurfaces: Geometry, Topology, and Algorithms"
			// is implemente here.
			// If a problematic C16 or C19 configuration shares the ambiguous face 
			// with another C16 or C19 configuration we simply invert the cube code
			// before looking up dual points. Doing this for these pairs ensures
			// manifold meshes.
			// But this removes the dualism to marching cubes.

			// check if we have a potentially problematic configuration
			uint8_t const direction = problematicConfigs[uint8_t(cubeCode)];
			// If the direction code is in {0,...,5} we have a C16 or C19 configuration.
			if (direction != 255) {
				// We have to check the neighboring cube, which shares the ambiguous
				// face. For this we decode the direction. This could also be done
				// with another lookup table.
				// copy current cube coordinates into an array.
				int32_t neighborCoords[] = { cx,cy,cz };
				// get the dimension of the non-zero coordinate axis
				unsigned int const component = direction >> 1;
				// get the sign of the direction
				int32_t delta = (direction & 1) == 1 ? 1 : -1;
				// modify the correspong cube coordinate
				neighborCoords[component] += delta;
				// have we left the volume in this direction?
				if (neighborCoords[component] >= 0 && neighborCoords[component] < (dims[component] - 1)) {
					// get the cube configuration of the relevant neighbor
					int neighborCubeCode = getCellCode(neighborCoords[0], neighborCoords[1], neighborCoords[2], iso);
					// Look up the neighbor configuration ambiguous face direction.
					// If the direction is valid we have a C16 or C19 neighbor.
					// As C16 and C19 have exactly one ambiguous face this face is
					// guaranteed to be shared for the pair.
					if (problematicConfigs[uint8_t(neighborCubeCode)] != 255) {
						// replace the cube configuration with its inverse.
						cubeCode ^= 0xff;
					}
				}
			}
		}
		for (int i = 0; i < 4; ++i)
			if (dualPointsList[cubeCode][i] & edge) {
				return dualPointsList[cubeCode][i];
			}
		return 0;
	}

	//------------------------------------------------------------------------------

	template<class T> inline
		void DualMC<T>::calculateDualPoint(int32_t const cx, int32_t const cy, int32_t const cz, VolumeDataType const iso, int const pointCode, Vertex& v) const {
		// initialize the point with lower voxel coordinates
		v.x = cx;
		v.y = cy;
		v.z = cz;

		// compute the dual point as the mean of the face vertices belonging to the
		// original marching cubes face
		Vertex p;
		p.x = 0;
		p.y = 0;
		p.z = 0;
		int points = 0;

		// sum edge intersection vertices using the point code
		if (pointCode & EDGE0) {
			p.x += ((float)iso - (float)data[gA(cx, cy, cz)]) / ((float)data[gA(cx + 1, cy, cz)] - (float)data[gA(cx, cy, cz)]);
			points++;
		}

		if (pointCode & EDGE1) {
			p.x += 1.0f;
			p.z += ((float)iso - (float)data[gA(cx + 1, cy, cz)]) / ((float)data[gA(cx + 1, cy, cz + 1)] - (float)data[gA(cx + 1, cy, cz)]);
			points++;
		}

		if (pointCode & EDGE2) {
			p.x += ((float)iso - (float)data[gA(cx, cy, cz + 1)]) / ((float)data[gA(cx + 1, cy, cz + 1)] - (float)data[gA(cx, cy, cz + 1)]);
			p.z += 1.0f;
			points++;
		}

		if (pointCode & EDGE3) {
			p.z += ((float)iso - (float)data[gA(cx, cy, cz)]) / ((float)data[gA(cx, cy, cz + 1)] - (float)data[gA(cx, cy, cz)]);
			points++;
		}

		if (pointCode & EDGE4) {
			p.x += ((float)iso - (float)data[gA(cx, cy + 1, cz)]) / ((float)data[gA(cx + 1, cy + 1, cz)] - (float)data[gA(cx, cy + 1, cz)]);
			p.y += 1.0f;
			points++;
		}

		if (pointCode & EDGE5) {
			p.x += 1.0f;
			p.z += ((float)iso - (float)data[gA(cx + 1, cy + 1, cz)]) / ((float)data[gA(cx + 1, cy + 1, cz + 1)] - (float)data[gA(cx + 1, cy + 1, cz)]);
			p.y += 1.0f;
			points++;
		}

		if (pointCode & EDGE6) {
			p.x += ((float)iso - (float)data[gA(cx, cy + 1, cz + 1)]) / ((float)data[gA(cx + 1, cy + 1, cz + 1)] - (float)data[gA(cx, cy + 1, cz + 1)]);
			p.z += 1.0f;
			p.y += 1.0f;
			points++;
		}

		if (pointCode & EDGE7) {
			p.z += ((float)iso - (float)data[gA(cx, cy + 1, cz)]) / ((float)data[gA(cx, cy + 1, cz + 1)] - (float)data[gA(cx, cy + 1, cz)]);
			p.y += 1.0f;
			points++;
		}

		if (pointCode & EDGE8) {
			p.y += ((float)iso - (float)data[gA(cx, cy, cz)]) / ((float)data[gA(cx, cy + 1, cz)] - (float)data[gA(cx, cy, cz)]);
			points++;
		}

		if (pointCode & EDGE9) {
			p.x += 1.0f;
			p.y += ((float)iso - (float)data[gA(cx + 1, cy, cz)]) / ((float)data[gA(cx + 1, cy + 1, cz)] - (float)data[gA(cx + 1, cy, cz)]);
			points++;
		}

		if (pointCode & EDGE10) {
			p.x += 1.0f;
			p.y += ((float)iso - (float)data[gA(cx + 1, cy, cz + 1)]) / ((float)data[gA(cx + 1, cy + 1, cz + 1)] - (float)data[gA(cx + 1, cy, cz + 1)]);
			p.z += 1.0f;
			points++;
		}

		if (pointCode & EDGE11) {
			p.z += 1.0f;
			p.y += ((float)iso - (float)data[gA(cx, cy, cz + 1)]) / ((float)data[gA(cx, cy + 1, cz + 1)] - (float)data[gA(cx, cy, cz + 1)]);
			points++;
		}

		// divide by number of accumulated points
		float invPoints = 1.0f / (float)points;
		p.x *= invPoints;
		p.y *= invPoints;
		p.z *= invPoints;

		// offset point by voxel coordinates
		v.x += p.x;
		v.y += p.y;
		v.z += p.z;
	}

	//------------------------------------------------------------------------------

	template<class T> inline
		QuadIndexType DualMC<T>::getSharedDualPointIndex(
			int32_t const cx, int32_t const cy, int32_t const cz,
			VolumeDataType const iso, DMCEdgeCode const edge,
			std::vector<Vertex>& vertices,
			PointToIndexMap& pointToIndex
		) {
		// create a key for the dual point from its linearized cell ID and point code
		DualPointKey key;
		key.linearizedCellID = gA(cx, cy, cz);
		key.pointCode = getDualPointCode(cx, cy, cz, iso, edge);

		// have we already computed the dual point?
		auto iterator = pointToIndex.find(key);
		if (iterator != pointToIndex.end()) {
			// just return the dual point index
			return iterator->second;
		}
		else {
			// create new vertex and vertex id
			QuadIndexType newVertexId = vertices.size();
			vertices.emplace_back();
			calculateDualPoint(cx, cy, cz, iso, key.pointCode, vertices.back());
			// insert vertex ID into map and also return it
			pointToIndex[key] = newVertexId;
			return newVertexId;
		}
	}

	//------------------------------------------------------------------------------

	template<class T> inline
		void DualMC<T>::build(
			VolumeDataType const* data,
			int32_t const dimX, int32_t const dimY, int32_t const dimZ,
			VolumeDataType const iso,
			bool const generateManifold,
			bool const generateSoup,
			std::vector<Vertex>& vertices,
			std::vector<Quad>& quads
		) {

		// set members
		this->dims[0] = dimX;
		this->dims[1] = dimY;
		this->dims[2] = dimZ;
		this->data = data;
		this->generateManifold = generateManifold;

		// clear vertices and quad indices
		vertices.clear();
		quads.clear();

		// generate quad soup or shared vertices quad list
		//if (generateSoup) {
			buildQuadSoup(iso, vertices, quads);
		//}
		///else {
			//buildSharedVerticesQuads(iso, vertices, quads);
		//}
	}

	//------------------------------------------------------------------------------

	template<class T> inline
		void DualMC<T>::buildQuadSoup(
			VolumeDataType const iso,
			std::vector<Vertex>& vertices,
			std::vector<Quad>& quads
		) {

		int32_t const reducedX = dims[0] - 2;
		int32_t const reducedY = dims[1] - 2;
		int32_t const reducedZ = dims[2] - 2;

		if (callback != NULL) callback(0);
		size_t total_progress = reducedZ;
		size_t progress = 0;

		// iterate voxels
#ifndef USE_PARALLEL
		for (int32_t z = 0; z < reducedZ; ++z) {
#define local_vertices vertices
#else
		static MutexType writeMutex;
		static MutexType callbackMutex;
		static tbb::affinity_partitioner ap;
		tbb::parallel_for(
			tbb::blocked_range<int32_t>(0, reducedZ),
			[&](const tbb::blocked_range<int32_t> r) {
				std::vector<Vertex> local_vertices;
				Vertex vertex0;
				Vertex vertex1;
				Vertex vertex2;
				Vertex vertex3;
				int pointCode;
				for (int32_t z = r.begin(); z < r.end(); z++) {
#endif
					for (int32_t y = 0; y < reducedY; ++y)
						for (int32_t x = 0; x < reducedX; ++x) {
							// construct quad for x edge
							if (z > 0 && y > 0) {
								// is edge intersected?
								bool const entering = data[gA(x, y, z)] < iso && data[gA(x + 1, y, z)] >= iso;
								bool const exiting = data[gA(x, y, z)] >= iso && data[gA(x + 1, y, z)] < iso;
								if (entering || exiting) {
									// generate quad
									pointCode = getDualPointCode(x, y, z, iso, EDGE0);
									calculateDualPoint(x, y, z, iso, pointCode, vertex0);

									pointCode = getDualPointCode(x, y, z - 1, iso, EDGE2);
									calculateDualPoint(x, y, z - 1, iso, pointCode, vertex1);

									pointCode = getDualPointCode(x, y - 1, z - 1, iso, EDGE6);
									calculateDualPoint(x, y - 1, z - 1, iso, pointCode, vertex2);

									pointCode = getDualPointCode(x, y - 1, z, iso, EDGE4);
									calculateDualPoint(x, y - 1, z, iso, pointCode, vertex3);

									if (entering) {
										local_vertices.emplace_back(vertex0);
										local_vertices.emplace_back(vertex1);
										local_vertices.emplace_back(vertex2);
										local_vertices.emplace_back(vertex3);
									}
									else {
										local_vertices.emplace_back(vertex0);
										local_vertices.emplace_back(vertex3);
										local_vertices.emplace_back(vertex2);
										local_vertices.emplace_back(vertex1);
									}
								}
							}

							// construct quad for y edge
							if (z > 0 && x > 0) {
								// is edge intersected?
								bool const entering = data[gA(x, y, z)] < iso && data[gA(x, y + 1, z)] >= iso;
								bool const exiting = data[gA(x, y, z)] >= iso && data[gA(x, y + 1, z)] < iso;
								if (entering || exiting) {
									// generate quad
									pointCode = getDualPointCode(x, y, z, iso, EDGE8);
									calculateDualPoint(x, y, z, iso, pointCode, vertex0);

									pointCode = getDualPointCode(x, y, z - 1, iso, EDGE11);
									calculateDualPoint(x, y, z - 1, iso, pointCode, vertex1);

									pointCode = getDualPointCode(x - 1, y, z - 1, iso, EDGE10);
									calculateDualPoint(x - 1, y, z - 1, iso, pointCode, vertex2);

									pointCode = getDualPointCode(x - 1, y, z, iso, EDGE9);
									calculateDualPoint(x - 1, y, z, iso, pointCode, vertex3);

									if (exiting) {
										local_vertices.emplace_back(vertex0);
										local_vertices.emplace_back(vertex1);
										local_vertices.emplace_back(vertex2);
										local_vertices.emplace_back(vertex3);
									}
									else {
										local_vertices.emplace_back(vertex0);
										local_vertices.emplace_back(vertex3);
										local_vertices.emplace_back(vertex2);
										local_vertices.emplace_back(vertex1);
									}
								}
							}

							// construct quad for z edge
							if (x > 0 && y > 0) {
								// is edge intersected?
								bool const entering = data[gA(x, y, z)] < iso && data[gA(x, y, z + 1)] >= iso;
								bool const exiting = data[gA(x, y, z)] >= iso && data[gA(x, y, z + 1)] < iso;
								if (entering || exiting) {
									// generate quad
									pointCode = getDualPointCode(x, y, z, iso, EDGE3);
									calculateDualPoint(x, y, z, iso, pointCode, vertex0);

									pointCode = getDualPointCode(x - 1, y, z, iso, EDGE1);
									calculateDualPoint(x - 1, y, z, iso, pointCode, vertex1);

									pointCode = getDualPointCode(x - 1, y - 1, z, iso, EDGE5);
									calculateDualPoint(x - 1, y - 1, z, iso, pointCode, vertex2);

									pointCode = getDualPointCode(x, y - 1, z, iso, EDGE7);
									calculateDualPoint(x, y - 1, z, iso, pointCode, vertex3);

									if (exiting) {
										local_vertices.emplace_back(vertex0);
										local_vertices.emplace_back(vertex1);
										local_vertices.emplace_back(vertex2);
										local_vertices.emplace_back(vertex3);
									}
									else {
										local_vertices.emplace_back(vertex0);
										local_vertices.emplace_back(vertex3);
										local_vertices.emplace_back(vertex2);
										local_vertices.emplace_back(vertex1);
									}
								}
							}
							//++progress;
						}
				}
#ifdef USE_PARALLEL
				{
					MutexType::scoped_lock lock(writeMutex);
					vertices.insert(vertices.end(), local_vertices.begin(), local_vertices.end());
					progress += (r.end() - r.begin());
				}
				if (callback != NULL) {
					int p = progress * 100 / (total_progress);
					CallbackMutexType::scoped_lock lock(callbackMutex);
					callback(p >= 100 ? 99 : p);
				}
			}, ap);
#else
#undef local_vertices
				progress += (r.end() - r.begin());
				callback(progress * 100 / total_progress);
#endif
				// generate triangle soup quads
				size_t const numQuads = vertices.size() / 4;
				quads.reserve(numQuads);
				for (size_t i = 0; i < numQuads; ++i) {
					quads.emplace_back(i * 4, i * 4 + 1, i * 4 + 2, i * 4 + 3);
				}
				if (callback != NULL) {
					callback(100);
				}
		}

			/*/------------------------------------------------------------------------------

			template<class T> inline
				void DualMC<T>::buildSharedVerticesQuads(
					VolumeDataType const iso,
					std::vector<Vertex>& vertices,
					std::vector<Quad>& quads
				) {


				int32_t const reducedX = dims[0] - 2;
				int32_t const reducedY = dims[1] - 2;
				int32_t const reducedZ = dims[2] - 2;

				QuadIndexType i0, i1, i2, i3;

				//pointToIndex.clear();

				//ProgressBar progress(reducedX * reducedY * reducedZ, PROGRESS_BAR_COLUMN);

				// iterate voxels
		#ifndef USE_PARALLEL
				for (int32_t z = 0; z < reducedZ; ++z) {
		//#define local_quads quads
		//#define local_vertices vertices
		#else
				tbb::spin_mutex writeMutex;
				static tbb::affinity_partitioner ap;
				tbb::parallel_for(
				tbb::blocked_range<int32_t>(0, reducedZ),
				[&](const tbb::blocked_range<int32_t> r) {
					std::vector<Vertex> local_vertices;
					std::vector<Quad> local_quads;
					PointToIndexMap pointToIndex;
					for (int32_t z = r.begin(); z < r.end(); z++) {
		#endif
						for (int32_t y = 0; y < reducedY; ++y)
							for (int32_t x = 0; x < reducedX; ++x) {
								// construct quads for x edge
								if (z > 0 && y > 0) {
									bool const entering = data[gA(x, y, z)] < iso && data[gA(x + 1, y, z)] >= iso;
									bool const exiting = data[gA(x, y, z)] >= iso && data[gA(x + 1, y, z)] < iso;
									if (entering || exiting) {
										// generate quad
										i0 = getSharedDualPointIndex(x, y, z, iso, EDGE0, local_vertices, pointToIndex);
										i1 = getSharedDualPointIndex(x, y, z - 1, iso, EDGE2, local_vertices, pointToIndex);
										i2 = getSharedDualPointIndex(x, y - 1, z - 1, iso, EDGE6, local_vertices, pointToIndex);
										i3 = getSharedDualPointIndex(x, y - 1, z, iso, EDGE4, local_vertices, pointToIndex);

										if (entering) {
											local_quads.emplace_back(i0, i1, i2, i3);
										}
										else {
											local_quads.emplace_back(i0, i3, i2, i1);
										}
									}
								}

								// construct quads for y edge
								if (z > 0 && x > 0) {
									bool const entering = data[gA(x, y, z)] < iso && data[gA(x, y + 1, z)] >= iso;
									bool const exiting = data[gA(x, y, z)] >= iso && data[gA(x, y + 1, z)] < iso;
									if (entering || exiting) {
										// generate quad
										i0 = getSharedDualPointIndex(x, y, z, iso, EDGE8, local_vertices, pointToIndex);
										i1 = getSharedDualPointIndex(x, y, z - 1, iso, EDGE11, local_vertices, pointToIndex);
										i2 = getSharedDualPointIndex(x - 1, y, z - 1, iso, EDGE10, local_vertices, pointToIndex);
										i3 = getSharedDualPointIndex(x - 1, y, z, iso, EDGE9, local_vertices, pointToIndex);

										if (exiting) {
											local_quads.emplace_back(i0, i1, i2, i3);
										}
										else {
											local_quads.emplace_back(i0, i3, i2, i1);
										}
									}
								}

								// construct quads for z edge
								if (x > 0 && y > 0) {
									bool const entering = data[gA(x, y, z)] < iso && data[gA(x, y, z + 1)] >= iso;
									bool const exiting = data[gA(x, y, z)] >= iso && data[gA(x, y, z + 1)] < iso;
									if (entering || exiting) {
										// generate quad
										i0 = getSharedDualPointIndex(x, y, z, iso, EDGE3, local_vertices, pointToIndex);
										i1 = getSharedDualPointIndex(x - 1, y, z, iso, EDGE1, local_vertices, pointToIndex);
										i2 = getSharedDualPointIndex(x - 1, y - 1, z, iso, EDGE5, local_vertices, pointToIndex);
										i3 = getSharedDualPointIndex(x, y - 1, z, iso, EDGE7, local_vertices, pointToIndex);

										if (exiting) {
											local_quads.emplace_back(i0, i1, i2, i3);
										}
										else {
											local_quads.emplace_back(i0, i3, i2, i1);
										}
									}
								}
								// progress bar
								//++progress;
							}
						//progress.display();
					}
		#ifdef USE_PARALLEL
					{
						tbb::spin_mutex::scoped_lock lock(writeMutex);
						QuadIndexType last_index = vertices.size();
						for (Quad& q : local_quads) {
							q.i0 += last_index;
							q.i1 += last_index;
							q.i2 += last_index;
							q.i3 += last_index;
						}
						vertices.insert(vertices.end(), local_vertices.begin(), local_vertices.end());
						quads.insert(quads.end(), local_quads.begin(), local_quads.end());
					}
				}, ap);
				std::cout << vertices.size() << "," << quads.size() << std::endl;
		#endif
				//progress.done();
			}*/

#include "dualmc_tables.tpp"

	} // END: namespace dualmc
#endif // DUALMC_H_INCLUDED