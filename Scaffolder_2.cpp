#include <igl/writePLY.h>
#include <igl/readSTL.h>
#include <igl/copyleft/marching_cubes.h>
#include <igl/fast_winding_number.h>
#include <Eigen/Core>
#include "Scaffolder_2.h"

int main(int argc, char* argv[])
{
    // Define default parameters
    bool verbose = true;
    bool use_smooth_shape = false;
    bool use_smooth_mesh = false;
    bool is_test = false;
    int precision = 6;
    uint16_t border_offset = 2;
    double thickness = 0.5;
    size_t grid_size = 100;
    double coff = pi;
    // file parameters
    std::string filename = "out";
    std::string format = "ply";

    // Inline mesh of a cube
    Eigen::MatrixXd V1,V2;
    Eigen::MatrixXi F1,F2;
    {
        Eigen::MatrixXd N;
        igl::readSTL("2.stl", V1, F1, N);
    }
    
    const Eigen::RowVector3d V1min = V1.colwise().minCoeff();
    const Eigen::RowVector3d V1max = V1.colwise().maxCoeff();
    const Eigen::RowVector3d L = V1max - V1min;
    const Eigen::RowVector3d delta = L / grid_size;
    Iso_cuboid_condition condition(V1min(0), V1min(1), V1min(2), L(0), L(1), L(2));
    Implicit_function fn(condition, schwarzp, coff, thickness);
    const Eigen::RowVector3d V1min1 = V1min - border_offset * delta;
    grid_size += 2*(size_t)(border_offset);

    if (verbose) std::cout << "Bounding Box: " << V1min << ' ' << V1max << std::endl;
    if (verbose) std::cout << "Length: " << L << std::endl;
    //if (verbose) std::cout << "delta: " << delta << std::endl;
    if (verbose) std::cout << "[Generating grid] " << grid_size * grid_size * grid_size << ' (' << grid_size << '*' << grid_size << '*' << grid_size << ') ';
    Eigen::MatrixXd GV(grid_size * grid_size * grid_size, 3);
    Eigen::VectorXd Fxyz(grid_size * grid_size * grid_size);
    for (size_t k = 0; k < grid_size; k++) {
        const double z = V1min1(2) + k * delta(2);
        for (size_t j = 0; j < grid_size; j++) {
            const double y = V1min1(1) + j * delta(1);
            for (size_t i = 0; i < grid_size; i++) {
                const double x = V1min1(0) + i * delta(0);
                const size_t index = i + grid_size * (j + grid_size * k);
                GV.row(index) = Eigen::RowVector3d(x, y, z);
            }
        }
    }
    if (verbose) std::cout << "OK" << std::endl;

    if (verbose) std::cout << "[Calculating Wind number] ";
    Eigen::VectorXd W;
    {
        igl::FastWindingNumberBVH bvh;
        igl::fast_winding_number(V1, F1, 2, bvh);
        igl::fast_winding_number(bvh, 2, GV, W);
    }
    if (verbose) std::cout << "OK" << std::endl;
    if (verbose) std::cout << "-- Info: W>=0.5 " << (W.array() >= 0.5).count() << std::endl;

    if (verbose) std::cout << "[Generating isosurface Fxyz] ";
    for (size_t k = 0; k < grid_size; k++) {
        const double z = V1min1(2) + k * delta(2);
        for (size_t j = 0; j < grid_size; j++) {
            const double y = V1min1(1) + j * delta(1);
            for (size_t i = 0; i < grid_size; i++) {
                const double x = V1min1(0) + i * delta(0);
                const size_t index = i + grid_size * (j + grid_size * k);
                Fxyz(index) = W(index) >= 0.5 ? fn(x, y, z) : 1.0;
            }
        }
    }
    if (verbose) std::cout << "OK" << std::endl;
    
    if (verbose) std::cout << "[Marching Cube] ";
    igl::copyleft::marching_cubes(Fxyz, GV, grid_size, grid_size, grid_size, V2, F2);
    if (verbose) std::cout << "OK" << std::endl;

    
    if (verbose) std::cout << "[Writing PLY] ";
    igl::writePLY("out.ply", V2, F2);
    if (verbose) std::cout << "OK" << std::endl;
    /*
    if (verbose) std::cout << "[Writing OBJ] ";
    igl::writeOBJ("out.obj", SV, SF);
    if (verbose) std::cout << "OK" << std::endl;
    */
    return 0;
}