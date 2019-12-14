#include <igl/writePLY.h>
#include <igl/readSTL.h>
#include <igl/copyleft/marching_cubes.h>
#include <igl/fast_winding_number.h>
#include <Eigen/Core>
#include "Implicit_function.h"
#include "argparse.h"
#include "dualmc/dualmc.h"
#include "Mesh.h"

int main(int argc, const char* argv[])
{
    ArgumentParser parser("Scaffolder v1.0 - generate 3D scaffold from STL file");
    parser.add_argument("-q", "--quiet", "Disable verbose output", false);
    parser.add_argument("-f", "--format", "Output format (OFF,PLY,STL,OBJ) [default: ply]", false);
    parser.add_argument("-i", "--input", "Input file (STL)", true);
    parser.add_argument("-o", "--output", "Output filename without extension [default: out]", false);
    parser.add_argument("-g", "--grid", "Grid size [default: 100]", false);
    parser.add_argument("--thickness", "Thickness [default: 1.0]", false);
    parser.add_argument("--border_offset", "default:2", false);
    parser.add_argument("--coff", "default:4*PI", false);
    parser.add_argument("--minimum_diameter", "used for removing small orphaned (between 0-1) [default: 0.25]", false);
    parser.add_argument("--surface", "default: schwarzp", false);
    try {
        parser.parse(argc, argv);
    }
    catch (const ArgumentParser::ArgumentNotFound & ex) {
        std::cout << ex.what() << std::endl;
        return 0;
    }
    if (parser.is_help()) return 0;
    // Define default parameters
    bool verbose = true;
    bool use_smooth_mesh = false;
    uint16_t border_offset = 2;
    double thickness = 0.5;
    size_t grid_size = 100;
    double coff = pi;
    double minimum_diameter = 0.25;
    Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
    // file parameters
    std::string filename = "out";
    std::string format = "ply";
    std::string input_file = parser.get<std::string>("input");
    std::string surface = "schwarzp";

    // get optional parameters
    if (parser.exists("quiet")) verbose = parser.get<bool>("quiet");
    if (parser.exists("file")) filename = parser.get<std::string>("filename");
    if (parser.exists("format")) format = parser.get<std::string>("format");
    if (parser.exists("thickness")) thickness = parser.get<double>("thickness");
    if (parser.exists("grid_size")) grid_size = parser.get<size_t>("grid_size");
    if (parser.exists("border_offset")) border_offset = parser.get<uint16_t>("border_offset");
    if (parser.exists("coff")) coff = parser.get<double>("coff");
    if (parser.exists("minimum_diameter")) minimum_diameter = parser.get<double>("minimum_diameter");
    if (parser.exists("surface")) surface = parser.get<std::string>("surface");

    // lowercase format
    to_lower(format);
    if (format != "ply" && format != "obj" && format != "stl" && format != "off") {
        std::cout << "Invalid format: " << format << std::endl;
        return -1;
    }

    // Print parameters
    if (verbose) {
        std::cout << "[Scaffolder v1.0]" << std::endl
            << "-- Input file: " << input_file << std::endl
            << "-- Output file: " << filename << '.' << format << std::endl
            << "-- Coff: " << coff << std::endl
            << "-- Thickness: " << thickness << std::endl
            << "-- Grid size: " << grid_size << std::endl
            << "-- Border Offset: " << border_offset << std::endl
            << "-- Minimum diameter: " << 100 * minimum_diameter << "%" << std::endl;
    }

    // Stage 1:
    TMesh mesh;
    {
        Eigen::MatrixXd V, Fxyz;
        Eigen::MatrixXi F;
        Eigen::RowVector3d V1min1, delta;
        {
            // Read FROM STL
            Eigen::MatrixXd V1;
            Eigen::MatrixXi F1;
            {
                Eigen::MatrixXd N;
                igl::readSTL(input_file, V1, F1, N);
            }
            // Calculate the grid size parameters
            const Eigen::RowVector3d V1min = V1.colwise().minCoeff();
            const Eigen::RowVector3d V1max = V1.colwise().maxCoeff();
            const Eigen::RowVector3d L = V1max - V1min;
            delta = L / grid_size;
            // Create border offset from the original box
            V1min1 = V1min - border_offset * delta;
            grid_size += 2 * (size_t)(border_offset);
            // Create iso cuboid condition with boundary box
            //Iso_cuboid_condition condition(V1min(0), V1min(1), V1min(2), L(0), L(1), L(2));
            Function_3& isosurface = get_surface_function(surface);
            Implicit_function fn(isosurface, coff, thickness);
            if (verbose) std::cout << "Bounding Box: " << V1min.format(CleanFmt) << ' ' << V1max.format(CleanFmt) << std::endl;
            if (verbose) std::cout << "Length: " << L << std::endl;
            if (verbose) std::cout << "[Generating grid] " << grid_size * grid_size * grid_size << " (" << grid_size << '*' << grid_size << '*' << grid_size << ") ";
            Eigen::MatrixXd GV(grid_size * grid_size * grid_size, 3);
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
            //if (verbose) std::cout << "-- Info: W>=0.5 " << (W.array() >= 0.5).count() << std::endl;

            if (verbose) std::cout << "[Generating isosurface Fxyz] ";
            Fxyz.resize(grid_size * grid_size * grid_size, 1);
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
        }

        if (verbose) std::cout << "[Marching Cube] ";
        {
            dualmc::DualMC<double> builder;
            std::vector<dualmc::Vertex> mc_vertices;
            std::vector<dualmc::Quad> mc_quads;
            builder.build((double const*)Fxyz.data(), grid_size, grid_size, grid_size, 0, true, false, mc_vertices, mc_quads);
            TMesh::VertexIterator vi = vcg::tri::Allocator<TMesh>::AddVertices(mesh, mc_vertices.size());
            TMesh::FaceIterator fi = vcg::tri::Allocator<TMesh>::AddFaces(mesh, mc_quads.size()*2);
            std::vector<TMesh::VertexPointer> vp(mc_vertices.size());
            for (size_t i = 0, len = mc_vertices.size(); i < len; i++, ++vi) {
                vp[i] = &(*vi);
                vi->P() = TMesh::CoordType(
                    V1min1(0) + mc_vertices[i].x * delta(0),
                    V1min1(1) + mc_vertices[i].y * delta(1),
                    V1min1(2) + mc_vertices[i].z * delta(2)
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
            vcg::tri::Clean<TMesh>::RemoveDuplicateFace(mesh);
            vcg::tri::Clean<TMesh>::RemoveDuplicateVertex(mesh);
            vcg::tri::UpdateTopology<TMesh>::FaceFace(mesh);

        }
        if (verbose) std::cout << "OK" << std::endl;
    }

    // Stage 2: Clearning
    {
        if (verbose) std::cout << "[libVCG Cleaning] ";
        vcg::tri::Clean<TMesh>::RemoveZeroAreaFace(mesh);
        vcg::tri::Clean<TMesh>::RemoveUnreferencedVertex(mesh);
        vcg::tri::UpdateTopology<TMesh>::FaceFace(mesh);
        vcg::tri::Clean<TMesh>::RemoveNonManifoldFace(mesh);
        vcg::tri::Clean<TMesh>::RemoveUnreferencedVertex(mesh);
        vcg::tri::UpdateTopology<TMesh>::FaceFace(mesh);
        vcg::tri::UpdateBounding<TMesh>::Box(mesh);
        vcg::tri::Clean<TMesh>::RemoveSmallConnectedComponentsDiameter(mesh, minimum_diameter*mesh.bbox.Diag());
        vcg::tri::Clean<TMesh>::RemoveUnreferencedVertex(mesh);
        vcg::tri::UpdateTopology<TMesh>::FaceFace(mesh);
        if (mesh.fn > 0) {
            vcg::tri::UpdateNormal<TMesh>::PerFaceNormalized(mesh);
            vcg::tri::UpdateNormal<TMesh>::PerVertexAngleWeighted(mesh);
        }
        if (verbose) {
            std::cout
                << "OK" << std::endl
                << "-- IsWaterTight: " << vcg::tri::Clean<TMesh>::IsWaterTight(mesh) << std::endl
                << "-- Minimum_diameter: " << minimum_diameter * mesh.bbox.Diag() << std::endl;
        }
        
        if (verbose) std::cout << "[Writing file] ";
        filename.append(".");
        filename.append(format);
        if (format == "ply") {
            vcg::tri::io::ExporterPLY<TMesh>::Save(mesh, filename.c_str());
        }
        else if (format == "obj") {
            vcg::tri::io::ExporterOBJ<TMesh>::Save(mesh, filename.c_str(), 0);
        }
        else if (format == "off") {
            vcg::tri::io::ExporterOFF<TMesh>::Save(mesh, filename.c_str(), 0);
        }
        else if (format == "stl") {
            vcg::tri::io::ExporterSTL<TMesh>::Save(mesh, filename.c_str(), 0);
        }
        if (verbose) std::cout << "OK" << std::endl << "[Finished]" << std::endl;
    }
    
    return 0;
}