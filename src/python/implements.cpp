#include "PyScaffolder.hpp"
#include "OptimalSlice.hpp"
#include "MeshOperation.h"
#include "implicit_function.h"
#include "QuadricSimplification.h"
#include "dualmc/dualmc.h"

#define LOG std::cout
#define IS_INSTANCE_ARRAYLIST(x) (py::isinstance<py::list>(x) || py::isinstance<py::tuple>(x) || py::isinstance<py::array>(x))

using namespace PyScaffolder;
using namespace optimal_slice;
namespace py = pybind11;

// format Eigen object to readable text
Eigen::IOFormat CleanFmt(4, Eigen::DontAlignCols, ", ", "\n", "[", "]");
Eigen::IOFormat CSVFmt(-1, Eigen::DontAlignCols, ", ", ", ");

void clean_mesh(TMesh &mesh) {
    // It's a good practice to clean STL mesh firstly (VCGLib)
    mesh.face.EnableFFAdjacency();
    vcg::tri::Clean<TMesh>::RemoveDuplicateFace(mesh);
    vcg::tri::Clean<TMesh>::RemoveDuplicateVertex(mesh);
    vcg::tri::Clean<TMesh>::RemoveUnreferencedVertex(mesh);
    vcg::tri::UpdateTopology<TMesh>::FaceFace(mesh);
    vcg::tri::Clean<TMesh>::RemoveZeroAreaFace(mesh);
    vcg::tri::Clean<TMesh>::RemoveNonManifoldFace(mesh);
    vcg::tri::Clean<TMesh>::RemoveUnreferencedVertex(mesh);
    vcg::tri::UpdateTopology<TMesh>::FaceFace(mesh);
    mesh.face.DisableFFAdjacency();
}

MeshInfo PyScaffolder::generate_scaffold(
	Eigen::MatrixXd V1,
	Eigen::MatrixXi F1,
	Parameter params,
	const std::function<void(int)>& callback
) {
	MeshInfo mesh_info;
    TMesh mesh;
	double volume1 = EPSIL, volume2 = EPSIL;
	double area1 = EPSIL, area2 = EPSIL;
	{
		Eigen::MatrixXd V, Fxyz, IFxyz;
		Eigen::MatrixXi F;
		Eigen::RowVector3d V1min1;
		Eigen::RowVector3i grid_size;
		
        Eigen::RowVector3d V1min, V1max, L, delta;
        double grid_delta;
        {
            TMesh stl;
            eigen_vector_to_mesh(V1, F1, stl);
            // It's a good practice to clean STL mesh firstly (VCGLib)
            stl.face.EnableFFAdjacency();
            vcg::tri::Clean<TMesh>::RemoveDuplicateFace(stl);
            vcg::tri::Clean<TMesh>::RemoveDuplicateVertex(stl);
            vcg::tri::Clean<TMesh>::RemoveUnreferencedVertex(stl);
            vcg::tri::UpdateTopology<TMesh>::FaceFace(stl);
            vcg::tri::Clean<TMesh>::RemoveZeroAreaFace(stl);
            vcg::tri::Clean<TMesh>::RemoveNonManifoldFace(stl);
            vcg::tri::Clean<TMesh>::RemoveUnreferencedVertex(stl);
            vcg::tri::UpdateTopology<TMesh>::FaceFace(stl);
            stl.face.DisableFFAdjacency();
            // Check that the mesh is neither point cloud (has only vertices) nor non-watertight (has some self-intersacting faces).
            int edgeNum = 0, edgeBorderNum = 0, edgeNonManifNum = 0;
            vcg::tri::Clean<TMesh>::CountEdgeNum(stl, edgeNum, edgeBorderNum, edgeNonManifNum);
            bool watertight = (edgeBorderNum == 0) && (edgeNonManifNum == 0);
            bool pointcloud = (stl.fn == 0 && stl.vn != 0);
            if (!pointcloud && watertight) {
                // Collect the volume and surface area of original mesh
                area1 = vcg::tri::Stat<TMesh>::ComputeMeshArea(stl);
                volume1 = vcg::tri::Stat<TMesh>::ComputeMeshVolume(stl);
            }
            else {
                std::stringstream error_message;
                error_message << "Error: Input file is not valid" << std::endl
                    << "-- Watertight: " << watertight << (watertight ? "[Valid]" : "[Invalid]") << std::endl
                    << "-- Point cloud: " << pointcloud << (!pointcloud ? "[Valid]" : "[Invalid]") << std::endl;
                throw std::runtime_error(error_message.str());
            }
            

            // Calculate the grid size parameters
            V1min = V1.colwise().minCoeff();
            V1max = V1.colwise().maxCoeff();
            L = V1max - V1min;
            delta = L / params.grid_size;
            grid_delta = delta.minCoeff();
            // Create border offset from the original box
            V1min1 = V1min - params.grid_offset * grid_delta * Eigen::RowVector3d::Ones();
            grid_size = (L / grid_delta).cast<int>() + 2 * params.grid_offset * Eigen::RowVector3i::Ones();
            if (params.verbose)
                LOG << "-- Bounding Box: "
                    << V1min.format(CleanFmt) << ' ' << V1max.format(CleanFmt) << std::endl
                    << "-- Length: " << L.format(CleanFmt) << std::endl
                    << "-- Grid delta: " << grid_delta << std::endl;
            // Convert input STL from VCG::TMesh to Eigen matrixXd V1
            mesh_to_eigen_vector(stl, V1, F1);
        }
		// Generate Grid index (x,y,z)
        if (params.verbose) LOG << "[Generating grid] ";
		Eigen::MatrixXd GV(grid_size.prod(), 3);
		for (size_t k = 0; k < grid_size(2); k++) {
			const double z = V1min1(2) + k * grid_delta;
			for (size_t j = 0; j < grid_size(1); j++) {
				const double y = V1min1(1) + j * grid_delta;
				for (size_t i = 0; i < grid_size(0); i++) {
					const double x = V1min1(0) + i * grid_delta;
					const size_t index = i + grid_size(0) * (j + grid_size(1) * k);
					GV.row(index) = Eigen::RowVector3d(x, y, z);
				}
			}
		}
        if (params.verbose)
            LOG << "OK" << std::endl
                << "-- Grid size: " << grid_size.prod() << " " << grid_size.format(CleanFmt) << std::endl;
        
        if (callback != NULL) callback(1);

        Eigen::VectorXd W, D, Signed_Distance;
        {
            igl::FastWindingNumberBVH bvh;
            igl::AABB<Eigen::MatrixXd, 3> tree;
            tree.init(V1, F1);
            igl::fast_winding_number(V1, F1, 2, bvh);
            igl::fast_winding_number(bvh, 2, GV, W);
            igl::signed_distance_fast_winding_number(GV, V1, F1, tree, bvh, Signed_Distance);
            Signed_Distance *= -1;
        }
        if (params.verbose)
        LOG << "-- Sign Distance: [ " << Signed_Distance.minCoeff() << ", " << Signed_Distance.maxCoeff() << "]  "
            << "Wind: [ " << W.minCoeff() << ", " << W.maxCoeff() << "]" << std::endl;
        
        if (callback != NULL) callback(10);
        
        sol::state lua;
        Function* surface;
        sol::function f;

        // Initialize surface
        bool isLuaFunction = false;
        if (util::PathGetExtension(params.surface_name) == ".lua") {
            lua.open_libraries(sol::lib::base, sol::lib::math, sol::lib::string);
            set_shorten_function(lua);
            sol::load_result lua_file = lua.load_file(params.surface_name);
            if (!lua_file.valid()) {
                sol::error err = lua_file;
                std::cerr << "[Lua Error] " << err.what() << std::endl;
                exit(-1);
            }
            else {
                lua["params"] = lua.create_table_with(
                    "coff", params.coff,
                    "isolevel", params.isolevel,
                    "k_slice", params.k_slice,
                    "k_polygon", params.k_polygon
                );
                lua["bbox"] = lua.create_table_with(
                    "length", sol::as_table(std::vector<double>{L(0), L(1), L(2)}),
                    "min", sol::as_table(std::vector<double>{V1min(0), V1min(1), V1min(2)}),
                    "max", sol::as_table(std::vector<double>{V1max(0), V1max(1), V1max(2)}),
                    "grid_density", sol::as_table(std::vector<int>{grid_size(0), grid_size(1), grid_size(2)})
                );
                lua.set_function("winding", [W, V1min1, grid_delta, grid_size](double x, double y, double z) -> double
                    {
                        auto i = floor((x - V1min1(0)) / grid_delta);
                        auto j = floor((y - V1min1(1)) / grid_delta);
                        auto k = floor((z - V1min1(2)) / grid_delta);
                        const size_t index = i + grid_size(0) * (j + grid_size(1) * k);
                        if (index < 0 || index >= (size_t)grid_size(0) * grid_size(1) * grid_size(2))
                            return (double)-1;
                        return W(index);
                    });
                lua.set_function("signed_distance", [Signed_Distance, V1min1, grid_delta, grid_size](double x, double y, double z) -> double
                    {
                        auto i = floor((x - V1min1(0)) / grid_delta);
                        auto j = floor((y - V1min1(1)) / grid_delta);
                        auto k = floor((z - V1min1(2)) / grid_delta);
                        const size_t index = i + grid_size(0) * (j + grid_size(1) * k);
                        if (index < 0 || index >= (size_t)grid_size(0) * grid_size(1) * grid_size(2))
                            return (double)-1;
                        return Signed_Distance(index);
                    });
                lua_file();
                f = lua["surface"];
                surface = new LuaFunction(f);
                isLuaFunction = true;
            }
        }
        else {
            surface = isosurface(params.surface_name, params.isolevel);
        }

        Implicit_function fn(surface, params.coff);
        // test for fn function, in case of lua function
        // if it's invalid, the exception will occur
        fn(0.1, 0.23, 0.45);

        Fxyz.resize(grid_size.prod(), 1);
        IFxyz.resize(grid_size.prod(), 1);
        // voxelize the surface
        for (size_t k = 0; k < grid_size(2); k++) {
            const double z = V1min1(2) + k * grid_delta;
            for (size_t j = 0; j < grid_size(1); j++) {
                const double y = V1min1(1) + j * grid_delta;
                for (size_t i = 0; i < grid_size(0); i++) {
                    const double x = V1min1(0) + i * grid_delta;
                    const size_t index = i + grid_size(0) * (j + grid_size(1) * k);
                    const double s = Signed_Distance(index);
                    const double w = W(index);
                    const double fv = fn(x, y, z);
                    if (params.is_intersect && w < 0.8) {
                        Fxyz(index) = -1.0;
                        IFxyz(index) = -1.0;
                    }
                    else if (params.is_intersect && params.shell > 0 && s <= params.shell) {
                        Fxyz(index) = s;
                        IFxyz(index) = s;
                    }
                    else if (
                        !params.is_intersect &&
                        (
                            i < params.grid_offset ||
                            i >= grid_size(0) - params.grid_offset ||
                            j < params.grid_offset ||
                            j >= grid_size(1) - params.grid_offset ||
                            k < params.grid_offset ||
                            k >= grid_size(2) - params.grid_offset
                            )
                        ) {
                        Fxyz(index) = -1.0;
                        IFxyz(index) = -1.0;
                    }
                    else {
                        Fxyz(index) = fv;
                        IFxyz(index) = -fv;
                    }
                }
            }
        }
        if (callback != NULL) callback(20);
        // Create scaffolder and inverse mesh by dual marching cube
        if (params.is_build_inverse) {
            marching_cube(mesh, IFxyz, grid_size, V1min1, grid_delta, false, false);
        }
        else {
            marching_cube(mesh, Fxyz, grid_size, V1min1, grid_delta, false, false);
        }
        if (callback != NULL) callback(50);
        clean_mesh(mesh, params.minimum_diameter, params.smooth_step);
        if (params.fix_self_intersect > 0) fix_self_intersect_mesh(mesh, params.minimum_diameter, params.fix_self_intersect);

        if (callback != NULL) callback(60);
        if (params.qsim_percent > 0) {
            vcg::tri::mesh_quad_simplification(mesh, params.qsim_percent, null_callback);
        }
        if (callback != NULL) callback(90);

        bool is_manifold = is_mesh_manifold(mesh);
        if (!is_manifold) {
            is_manifold = fix_non_manifold_vertices(mesh, params.minimum_diameter, 5);
            is_manifold = is_manifold && fix_non_manifold_edges(mesh, params.minimum_diameter, 5);
        }
        if (!is_manifold) {
            if (params.verbose)
                LOG << "[Warning] Mesh is not two manifold" << std::endl;
        }
        else {
            area2 = vcg::tri::Stat<TMesh>::ComputeMeshArea(mesh);
            volume2 = vcg::tri::Stat<TMesh>::ComputeMeshVolume(mesh);
            mesh_info.porosity = 1 - abs(volume2 / volume1);
            mesh_info.surface_area = area2;
            mesh_info.surface_area_ratio = abs(area2 / area1);
        }
	}
    clean_mesh(mesh);
    mesh_to_eigen_vector(mesh, mesh_info.v, mesh_info.f);
    if (callback != NULL) callback(100);
	return mesh_info;
}

std::tuple<Eigen::MatrixXd, Eigen::MatrixXi> PyScaffolder::marching_cubes(
    Eigen::VectorXd& Fxyz,
    py::object& grid_size,
    py::object& _delta,
    const std::vector<double>& Vmin,
    bool clean,
    const std::function<void(int)>& callback
)   {
    try {
        dualmc::DualMC<double> builder;
        std::vector<dualmc::Vertex> mc_vertices;
        std::vector<dualmc::Quad> mc_quads;
        if (callback != NULL) builder.callback = callback;
        // Type conversion
        std::array<int32_t, 3> g;
        if (IS_INSTANCE_ARRAYLIST(grid_size)) {
            g = py::cast< std::array<int32_t, 3> >(grid_size);
        }
        else {
            int32_t gs = py::cast<int32_t>(grid_size);
            g[0] = gs;
            g[1] = gs;
            g[2] = gs;
        }
        std::array<double, 3> delta;
        if (IS_INSTANCE_ARRAYLIST(_delta)) {
            delta = py::cast< std::array<double, 3> >(_delta);
        }
        else {
            double d = py::cast<double>(_delta);
            delta[0] = d;
            delta[1] = d;
            delta[2] = d;
        }

        // Dual-Marching cubes
        builder.build((double const*)Fxyz.data(), g[0], g[1], g[2], 0, true, true, mc_vertices, mc_quads);

        Eigen::MatrixXd v;
        Eigen::MatrixXi f;
        if (clean) {
            TMesh mesh;
            TMesh::VertexIterator vi = vcg::tri::Allocator<TMesh>::AddVertices(mesh, mc_vertices.size());
            TMesh::FaceIterator fi = vcg::tri::Allocator<TMesh>::AddFaces(mesh, mc_quads.size() * 2);
            std::vector<TMesh::VertexPointer> vp(mc_vertices.size());
            for (size_t i = 0, len = mc_vertices.size(); i < len; i++, ++vi) {
                vp[i] = &(*vi);
                vi->P() = TMesh::CoordType(
                    Vmin[0] + mc_vertices[i].x * delta[0],
                    Vmin[1] + mc_vertices[i].y * delta[1],
                    Vmin[2] + mc_vertices[i].z * delta[2]
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
            clean_mesh(mesh);
            mesh_to_eigen_vector(mesh, v, f);
        }
        else {
            v.resize(mc_vertices.size(), 3);
            for (size_t i = 0, len = mc_vertices.size(); i < len; i++) {
                v.row(i) << Vmin[0] + mc_vertices[i].x * delta[0], Vmin[1] + mc_vertices[i].y * delta[1], Vmin[2] + mc_vertices[i].z * delta[2];
            }
            f.resize(mc_quads.size() * 2, 3);
            for (size_t i = 0, j = 0, len = mc_quads.size(); i < len; i++) {
                f.row(j) << mc_quads[i].i0, mc_quads[i].i1, mc_quads[i].i2;
                j++;
                f.row(j) << mc_quads[i].i2, mc_quads[i].i3, mc_quads[i].i0;
                j++;
            }
        }


        return make_tuple(v, f);
    }
    catch (std::exception& e) {
        if (callback != NULL) callback(100);
        throw std::runtime_error(e.what());
    }
}

PoreSize PyScaffolder::slice_test(Eigen::MatrixXd v, Eigen::MatrixXi f, size_t k_slice, size_t k_polygon, int direction, const std::function<void(int)>& callback) {
    PoreSize pore_size;
    TMesh mesh;
    eigen_vector_to_mesh(v, f, mesh);

    vcg::tri::UpdateBounding<TMesh>::Box(mesh);
    vcg::tri::Clean<TMesh>::RemoveDuplicateVertex(mesh);
    vcg::tri::Allocator<TMesh>::CompactEveryVector(mesh);
    vcg::tri::UpdateTopology<TMesh>::FaceFace(mesh);
    vcg::tri::Clean<TMesh>::RemoveDuplicateFace(mesh);
    vcg::tri::UpdateTopology<TMesh>::FaceFace(mesh);
    vcg::tri::Clean<TMesh>::RemoveZeroAreaFace(mesh);
    vcg::tri::UpdateTopology<TMesh>::FaceFace(mesh);
    vcg::tri::Clean<TMesh>::MergeCloseVertex(mesh, SLICE_PRECISION * 10);
    vcg::tri::UpdateTopology<TMesh>::FaceFace(mesh);
    if (callback != NULL) callback(20);

    std::vector<double> minFeret, maxFeret, shape[5];

    if (direction > 2) {
        for (int i = 0; i < 2; i++) {
            Slice s = incremental_slicing(mesh, k_slice, i);
            ContourSlice C = contour_construct(s, i);
            measure_feret_and_shape(C, 4, minFeret, maxFeret, shape);
            if (callback != NULL) callback(20 + 30 * (i + 1));
        }
    }
    else {
        Slice s = incremental_slicing(mesh, k_slice, direction);
        ContourSlice C = contour_construct(s, direction);
        measure_feret_and_shape(C, k_polygon, minFeret, maxFeret, shape);
    }

    if (minFeret.size() > 0) {
        std::sort(minFeret.begin(), minFeret.end());
        pore_size.minFeret = Eigen::VectorXd::Map(&minFeret[0], minFeret.size());
    }

    if (maxFeret.size() > 0) {
        std::sort(maxFeret.begin(), maxFeret.end());
        pore_size.maxFeret = Eigen::VectorXd::Map(&maxFeret[0], maxFeret.size());
    }
    if (callback != NULL) callback(100);
    return pore_size;
}