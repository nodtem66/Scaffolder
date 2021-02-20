#include "PyScaffolder.hpp"
#include "Scaffolder_2.h"

using namespace PyScaffolder;

MeshInfo PyScaffolder::generate_mesh(
	Eigen::MatrixXd V1,
	Eigen::MatrixXi F1,
	Parameter params,
	const std::function<void(int)>& callback
) {
	MeshInfo mesh_info;
	try {
        TMesh mesh;
		double volume1 = eps, volume2 = eps;
		double area1 = eps, area2 = eps;
		{
			Eigen::MatrixXd V, Fxyz, IFxyz;
			Eigen::MatrixXi F;
			Eigen::RowVector3d V1min1;
			Eigen::RowVector3i grid_size;
			double grid_delta;
			{
				// Calculate the grid size parameters
				const Eigen::RowVector3d V1min = V1.colwise().minCoeff();
				const Eigen::RowVector3d V1max = V1.colwise().maxCoeff();
				const Eigen::RowVector3d L = V1max - V1min;
				const Eigen::RowVector3d delta = L / params.grid_size;
				grid_delta = delta.minCoeff();
				// Create border offset from the original box
				V1min1 = V1min - params.grid_offset * grid_delta * Eigen::RowVector3d::Ones();
				grid_size = (L / grid_delta).cast<int>() + 2 * params.grid_offset * Eigen::RowVector3i::Ones();
			}
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
                    std::cout << "Error: Input file is not valid" << std::endl
                        << "-- Watertight: " << watertight << (watertight ? "[Valid]" : "[Invalid]") << std::endl
                        << "-- Point cloud: " << pointcloud << (!pointcloud ? "[Valid]" : "[Invalid]") << std::endl;
                    exit(-1);
                }
                // Convert input STL from VCG::TMesh to Eigen matrixXd V1
                mesh_to_eigen_vector(stl, V1, F1);
            }
			// Generate Grid index (x,y,z)
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
            if (callback != NULL) callback(1);

			Eigen::VectorXd W, D, Signed_Distance;
			{
				igl::FastWindingNumberBVH bvh;
				igl::AABB<Eigen::MatrixXd, 3> tree;
				tree.init(V1, F1);
				igl::fast_winding_number(V1, F1, 2, bvh);
				igl::fast_winding_number(bvh, 2, GV, W);
				igl::signed_distance_fast_winding_number(GV, V1, F1, tree, bvh, Signed_Distance);
			}
            if (callback != NULL) callback(10);

            sol::state lua;
            Function* surface;
            sol::function f;
            // Initialize surface
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
                    lua.set("w", params.coff);
                    lua.set("t", params.isolevel);
                    lua.set_function("winding", [W, V1min1, grid_delta, grid_size](double x, double y, double z) -> double
                        {
                            auto i = floor(x - V1min1(0) / grid_delta);
                            auto j = floor(y - V1min1(1) / grid_delta);
                            auto k = floor(z - V1min1(2) / grid_delta);
                            const size_t index = i + grid_size(0) * (j + grid_size(1) * k);
                            if (index < 0 || index >= (size_t)grid_size(0) * grid_size(1) * grid_size(2))
                                return (double)0;
                            return W(index);
                        });
                    lua.set_function("signed_distance", [Signed_Distance, V1min1, grid_delta, grid_size](double x, double y, double z) -> double
                        {
                            auto i = floor(x - V1min1(0) / grid_delta);
                            auto j = floor(y - V1min1(1) / grid_delta);
                            auto k = floor(z - V1min1(2) / grid_delta);
                            const size_t index = i + grid_size(0) * (j + grid_size(1) * k);
                            if (index < 0 || index >= (size_t)grid_size(0) * grid_size(1) * grid_size(2))
                                return (double)0;
                            return Signed_Distance(index);
                        });
                    lua_file();
                    f = lua["surface"];
                    surface = new LuaFunction(f);
                }
            }
            else {
                surface = isosurface(params.surface_name, params.isolevel);
            }
            Implicit_function fn(surface, params.coff);
            fn(0.1, 0.23, 0.45);
            Fxyz.resize(grid_size.prod(), 1);
            IFxyz.resize(grid_size.prod(), 1);
            for (size_t k = 0; k < grid_size(2); k++) {
                const double z = V1min1(2) + k * grid_delta;
                for (size_t j = 0; j < grid_size(1); j++) {
                    const double y = V1min1(1) + j * grid_delta;
                    for (size_t i = 0; i < grid_size(0); i++) {
                        const double x = V1min1(0) + i * grid_delta;
                        const size_t index = i + grid_size(0) * (j + grid_size(1) * k);
                        const double w = W(index);
                        // if the winding number < 0
                        // this region is classified as outer space
                        if (w < 0.8) {// Outside
                            Fxyz(index) = 1.0;
                            IFxyz(index) = 1.0;
                        }
                        // if the winding number > 1.1 
                        // we make this region as a solid wall
                        else if (w >= 1.1) {
                            Fxyz(index) = -1.0;
                            IFxyz(index) = -1.0;
                        }
                        else {// Inside
                            Fxyz(index) = fn(x, y, z);
                            IFxyz(index) = -Fxyz(index);
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
            fix_self_intersect_mesh(mesh, params.minimum_diameter, 5);

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
                std::cout << "[Warning] Mesh is not two manifold" << std::endl;
            }
            else {
                area2 = vcg::tri::Stat<TMesh>::ComputeMeshArea(mesh);
                volume2 = vcg::tri::Stat<TMesh>::ComputeMeshVolume(mesh);
                mesh_info.porosity = 1 - abs(volume2 / volume1);
                mesh_info.surface_area = area2;
                mesh_info.surface_area_ratio = abs(area2 / area1);
            }
		}
        mesh_to_eigen_vector(mesh, mesh_info.v, mesh_info.f);
        if (callback != NULL) callback(100);
	}
	catch (const std::exception& ex) {
		std::cout << "Exception: " << ex.what() << std::endl;
		std::cout << strerror(errno) << endl;
	}
	return mesh_info;
}

using namespace optimal_slice;


PoreSize PyScaffolder::slice_test(Eigen::MatrixXd v, Eigen::MatrixXi f, size_t k_slice, size_t k_polygon, int direction, const std::function<void(int)>& callback) {
    PoreSize pore_size;
    try {
        //std::cout << "VN: " << v.size() << " FN: " << f.size() << std::endl;
        //std::cout << "k_slice: " << k_slice << " direction: " << direction << std::endl;

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

    }
    catch (const std::exception& ex) {
        std::cout << "Exception: " << ex.what() << std::endl;
        std::cout << strerror(errno) << endl;
    }
    return pore_size;
}