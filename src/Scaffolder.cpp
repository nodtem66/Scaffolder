#include "Scaffolder.h"
#include "dualmc/dualmc.h"

ProgressBar qsim_progress(100, 40);
bool qsim_callback(int pos, const char* str) {
	if (pos >= 0 && pos <= 100) {
		qsim_progress.update(pos);
		qsim_progress.display();
	}
	if (pos >= 100)
		qsim_progress.done();
	return true;
}

ProgressBar marching_progress(100, 40);
void marching_callback(int pos) {
	if (pos >= 0 && pos <= 100) {
		marching_progress.update(pos);
		marching_progress.display();
	}
	if (pos >= 100)
		marching_progress.done();
}

// Flatten between 1D and 3D
// https://stackoverflow.com/questions/7367770/how-to-flatten-or-index-3d-array-in-1d-array
inline size_t indexFromIJK(size_t i, size_t j, size_t k, Eigen::RowVector3i grid_size) {
	return i + grid_size(0) * (j + grid_size(1) * k);
}

inline void indexToIJK(size_t index, Eigen::RowVector3i grid_size, index_type& r) {
	r.z = index / grid_size(0) / grid_size(1);
	index -= r.z * grid_size(0) * grid_size(1);
	r.y = index / grid_size(0);
	r.x = index % grid_size(0);
}

inline bool MarkAndSweepNeighbor(Eigen::VectorXd& W, index_type& index, Queue_t& queue, Eigen::RowVector3i grid_size, double value, bool findAbove) {
	bool isBorder = false;
	for (int8_t di = -1; di <= 1; di++) {
		for (int8_t dj = -1; dj <= 1; dj++) {
			for (int8_t dk = -1; dk <= 1; dk++) {
				if (di == 0 && dj == 0 && dk == 0) continue;
				const size_t id = indexFromIJK(index.x + di, index.y + dj, index.z + dk, grid_size);
				//std::cout << value << " " << W(id) << std::endl;
				if ((findAbove && W(id) >= value - EPSIL) || (!findAbove && W(id) <= value + EPSIL)) {
					isBorder = true;
					break;
				}
			}
			if (isBorder) break;
		}
		if (isBorder) break;
	}
	if (isBorder) {
		for (int8_t di = -1; di <= 1; di++) {
			for (int8_t dj = -1; dj <= 1; dj++) {
				for (int8_t dk = -1; dk <= 1; dk++) {
					if (di == 0 && dj == 0 && dk == 0) continue;
					const size_t id = indexFromIJK(index.x + di, index.y + dj, index.z + dk, grid_size);
					if (W(id) >= 0.5 && W(id) < 1.1) {
						queue.insert({ id, true });
					}
				}
			}
		}
	}
	return isBorder;
}

void marching_cube(TMesh& mesh, Eigen::MatrixXd& Fxyz, Eigen::RowVector3i grid_size, Eigen::RowVector3d& Vmin, double delta, bool verbose, bool dirty) {

	dualmc::DualMC<double> builder;
	std::vector<dualmc::Vertex> mc_vertices;
	std::vector<dualmc::Quad> mc_quads;
	if (verbose) {
		std::cout << "[Marching Cube] " << std::endl;
		builder.callback = marching_callback;
	}
	builder.build((double const*)Fxyz.data(), grid_size(0), grid_size(1), grid_size(2), 0, true, true, mc_vertices, mc_quads);
	TMesh::VertexIterator vi = vcg::tri::Allocator<TMesh>::AddVertices(mesh, mc_vertices.size());
	TMesh::FaceIterator fi = vcg::tri::Allocator<TMesh>::AddFaces(mesh, mc_quads.size() * 2);
	std::vector<TMesh::VertexPointer> vp(mc_vertices.size());
	for (size_t i = 0, len = mc_vertices.size(); i < len; i++, ++vi) {
		vp[i] = &(*vi);
		vi->P() = TMesh::CoordType(
			Vmin(0) + mc_vertices[i].x * delta,
			Vmin(1) + mc_vertices[i].y * delta,
			Vmin(2) + mc_vertices[i].z * delta
		);
	}
	for (size_t i = 0, len = mc_quads.size(); i < len; i++, ++fi) {
		fi->V(0) = vp[mc_quads[i].i0];
		fi->V(1) = vp[mc_quads[i].i1];
		fi->V(2) = vp[mc_quads[i].i2];
		++fi;
		fi->V(0) = vp[mc_quads[i].i2];
		fi->V(1) = vp[mc_quads[i].i3];
		fi->V(2) = vp[mc_quads[i].i0];
	}
	if (!dirty) {
		vcg::tri::Clean<TMesh>::RemoveDuplicateFace(mesh);
		vcg::tri::Clean<TMesh>::RemoveDuplicateVertex(mesh);
		vcg::tri::Clean<TMesh>::RemoveUnreferencedVertex(mesh);
	}

}

bool null_callback(int pos, const char* str) {
	return true;
}