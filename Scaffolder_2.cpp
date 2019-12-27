#include <igl/writePLY.h>
#include <igl/readSTL.h>
#include <igl/copyleft/marching_cubes.h>
#include <igl/fast_winding_number.h>
#include <Eigen/Core>
#include "Implicit_function.h"
#include "cxxopts.hpp"
#include "dualmc/dualmc.h"
#include "Mesh.h"

#define VERSION "v1.1"

typedef struct index_type {
    size_t x; size_t y; size_t z;
} index_type;
typedef std::map<size_t, bool> Queue_t;

inline size_t indexFromIJK(size_t i, size_t j, size_t k, size_t grid_size) {
    return i + grid_size * (j + grid_size * k);
}

inline void indexToIJK(size_t index, size_t grid_size, index_type& r) {
    r.z = index / (grid_size*grid_size);
    index -= r.z * grid_size * grid_size;
    r.y = index / grid_size;
    r.x = index % grid_size;
}

inline bool MarkAndSweepNeighbor(Eigen::VectorXd& W, index_type& index, Queue_t& queue, size_t grid_size, double value = 0.0, bool findAbove = false) {
    bool isBorder = false;
    for (int8_t di = -1; di <= 1; di++) {
        for (int8_t dj = -1; dj <= 1; dj++) {
            for (int8_t dk = -1; dk <= 1; dk++) {
                if (di == 0 && dj == 0 && dk == 0) continue;
                const size_t id = indexFromIJK(index.x + di, index.y + dj, index.z + dk, grid_size);
                //std::cout << value << " " << W(id) << std::endl;
                if ((findAbove && W(id) >= value - eps2) || (!findAbove && W(id) <= value + eps2)) {
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

int main(int argc, char* argv[])
{
    // Define default parameters
    bool verbose = true;
    bool dirty = false;
    uint16_t grid_offset = 2;
    uint16_t shell = 0;
    double thickness = 0.5;
    size_t grid_size = 100;
    double coff = pi;
    double minimum_diameter = 0.25;
    Eigen::IOFormat CleanFmt(4, Eigen::DontAlignCols, ", ", "\n", "[", "]");
    Eigen::IOFormat CSVFmt(-1, Eigen::DontAlignCols, ", ", ", ");
    // file parameters
    std::string filename = "out";
    std::string format = "ply";
    std::string surface = "schwarzp";
    std::string input_file = "";

    try {
        cxxopts::Options options("Scaffolder", "Scaffolder - generate 3D scaffold from STL file");
        options.positional_help("[option args]").show_positional_help();
        options.add_options()
            ("h,help", "Print help")
            ("q,quiet", "Disable verbose output")
            ("f,format", "Output format (OFF,PLY,STL,OBJ) [default: ply]", cxxopts::value<std::string>())
            ("i,input", "Input file (STL)", cxxopts::value<std::string>(), "FILE")
            ("o,output", "Output filename without extension [default: out]", cxxopts::value<std::string>(), "FILENAME")
            ("c,coff", "default:4*PI", cxxopts::value<double>(), "DOUBLE")
            ("s,shell", "[default:0]", cxxopts::value<uint16_t>(), "INT")
            ("n,surface", "rectlinear, schwarzp, schwarzd, gyroid, lidinoid, schoen_iwp, neovius, pwhybrid [default: schwarzp]", cxxopts::value<std::string>(), "NAME")
            ("t,thickness", "Thickness [default: 1.0]", cxxopts::value<double>(), "DOUBLE")
            ("g,grid_size", "Grid size [default: 100]", cxxopts::value<size_t>(), "INT")
            ("grid_offset", "[default:2]", cxxopts::value<uint16_t>(), "INT")
            ("dirty", "Disable autoclean")
            ("minimum_diameter", "used for removing small orphaned (between 0-1) [default: 0.25]", cxxopts::value<double>(), "DOUBLE");
        options.parse_positional({ "input", "output", "format" });
        bool isEmptyOption = (argc == 1);
        cxxopts::ParseResult result = options.parse(argc, argv);
        if (isEmptyOption || result.count("help")) {
            std::cout << options.help() << std::endl;
            return 0;
        }
        // Requirment
        if (result.count("input")) input_file = result["input"].as<std::string>();
        else {
            std::cout << "Missing Input file" << std::endl;
            return 1;
        }
        // get optional parameters
        if (result.count("quiet")) verbose = !result["quiet"].as<bool>();
        if (result.count("dirty")) dirty = result["dirty"].as<bool>();
        if (result.count("output")) filename = result["output"].as<std::string>();
        if (result.count("format")) format = result["format"].as<std::string>();
        if (result.count("thickness")) thickness = result["thickness"].as<double>();
        if (result.count("grid_size")) grid_size = result["grid_size"].as<size_t>();
        if (result.count("grid_offset")) grid_offset = result["grid_offset"].as<uint16_t>();
        if (result.count("coff")) coff = result["coff"].as<double>();
        if (result.count("minimum_diameter")) minimum_diameter = result["minimum_diameter"].as<double>();
        if (result.count("surface")) surface = result["surface"].as<std::string>();
        if (result.count("shell")) shell = result["shell"].as<uint16_t>();
        
        to_lower(surface);
        to_lower(format);
        
        if (format != "ply" && format != "obj" && format != "stl" && format != "off") {
            std::cout << "Invalid format: " << format << std::endl;
            return 1;
        }

        if (surface == "rectlinear") {
            thickness = 0;
        }
    }
    catch (const cxxopts::OptionException & ex) {
        std::cout << "Error parsing options: " << ex.what() << std::endl;
        return 1;
    }

    // Print parameters
    if (verbose) {
        std::cout << "[Scaffolder " << VERSION << "]" << std::endl
            << "-- Input file: " << input_file << std::endl
            << "-- Output file: " << filename << '.' << format << std::endl
            << "-- Surface: " << surface << std::endl
            << "-- Coff: " << coff << std::endl
            << "-- shell: " << shell << std::endl
            << "-- Thickness: " << thickness << std::endl
            << "-- Grid size: " << grid_size << std::endl
            << "-- Grid offset: " << grid_offset << std::endl
            << "-- Autoclean: " << (dirty ? "False" : "True") << std::endl
            << "-- Minimum diameter: " << 100 * minimum_diameter << "%" << std::endl;
    }

    // Stage 1:
    TMesh mesh;
    {
        Eigen::MatrixXd V, Fxyz;
        Eigen::MatrixXi F;
        Eigen::RowVector3d V1min1, delta;
        try {
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
            V1min1 = V1min - grid_offset * delta;
            grid_size += 2 * (size_t)(grid_offset);
            Implicit_function fn(isosurface(surface, coff), coff, thickness);
            if (verbose) std::cout << "-- Bounding Box: " << V1min.format(CleanFmt) << ' ' << V1max.format(CleanFmt) << std::endl;
            if (verbose) std::cout << "-- Length: " << L << std::endl;
            if (verbose) std::cout << "[Generating grid] ";
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
            if (verbose) 
                std::cout
                    << "OK" << std::endl
                    << "-- Grid size: " << grid_size * grid_size * grid_size << " (" << grid_size << 'x' << grid_size << 'x' << grid_size << ") " << std::endl;

            if (verbose) std::cout << "[Calculating Winding number] ";
            Eigen::VectorXd W;
            {
                igl::FastWindingNumberBVH bvh;
                igl::fast_winding_number(V1, F1, 2, bvh);
                igl::fast_winding_number(bvh, 2, GV, W);
            }
            if (verbose) std::cout << "OK" << std::endl;

            if (shell > 0) {
                if (verbose) std::cout << "[Generate Shell] ";
                Queue_t queue;
                uint16_t shell_index = shell + 1;
                for (size_t index = 0; index < grid_size * grid_size * grid_size; index++) {
                    if (W(index) >= 0.8) {
                        index_type id;
                        indexToIJK(index, grid_size, id);
                        if (MarkAndSweepNeighbor(W, id, queue, grid_size, 0.5)) {
                            W(index) = shell_index;
                        }
                    }
                }
                //if (verbose) std::cout << "-- [" << shell_index << "] " << queue.size() << std::endl;
                shell_index--;
                for (; shell_index > 1; shell_index--) {
                    Queue_t q(queue);
                    queue.clear();
                    //if (verbose) std::cout << "-- [" << shell_index << "] " << q.size() << std::endl;
                    for (Queue_t::iterator it = q.begin(); it != q.end(); ++it) {
                        // every neighbor
                        index_type id;
                        indexToIJK(it->first, grid_size, id);
                        if (MarkAndSweepNeighbor(W, id, queue, grid_size, shell_index, true)) {
                            W(it->first) = shell_index;
                        }
                    }
                }
                if (verbose) std::cout << "OK" << std::endl;
            }
            
            if (verbose) std::cout << "[Generating isosurface Fxyz] ";
            Fxyz.resize(grid_size * grid_size * grid_size, 1);
            for (size_t k = 0; k < grid_size; k++) {
                const double z = V1min1(2) + k * delta(2);
                for (size_t j = 0; j < grid_size; j++) {
                    const double y = V1min1(1) + j * delta(1);
                    for (size_t i = 0; i < grid_size; i++) {
                        const double x = V1min1(0) + i * delta(0);
                        const size_t index = i + grid_size * (j + grid_size * k);
                        const double w = W(index);
                        // Winding field
                        // Fxyz(index) = W(index);
                        if (w < 0.8) // Outside
                            Fxyz(index) = 1.0;
                        else if (w >= 1.1) // Shell
                            Fxyz(index) = -1;
                        else  // Inside
                            Fxyz(index) = fn(x, y, z);
                    }
                }
            }
            if (verbose) std::cout << "OK" << std::endl;
        }
        catch (const std::exception& ex) {
            std::cout << "Exception: " << ex.what() << std::endl;
            return 1;
        }

        {
            if (verbose) std::cout << "[Marching Cube] ";
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
            if (!dirty) {
                vcg::tri::Clean<TMesh>::RemoveDuplicateFace(mesh);
                vcg::tri::Clean<TMesh>::RemoveDuplicateVertex(mesh);
                vcg::tri::Clean<TMesh>::RemoveUnreferencedVertex(mesh);
                vcg::tri::UpdateTopology<TMesh>::FaceFace(mesh);
            }

            if (verbose) std::cout << "OK" << std::endl;
            if (verbose) std::cout << "-- Info: " << mc_vertices.size() << " vertices " << mc_quads.size() << " faces" << std::endl;
        }
    }

    // Stage 2: Cleaning
    {
        if (!dirty) {
            if (verbose) std::cout << "[libVCG Cleaning] ";
            vcg::tri::Clean<TMesh>::RemoveZeroAreaFace(mesh);
            vcg::tri::Clean<TMesh>::RemoveNonManifoldFace(mesh);
            vcg::tri::Clean<TMesh>::RemoveUnreferencedVertex(mesh);
            vcg::tri::UpdateTopology<TMesh>::FaceFace(mesh);
            vcg::tri::UpdateBounding<TMesh>::Box(mesh);
            vcg::tri::Clean<TMesh>::RemoveSmallConnectedComponentsDiameter(mesh, minimum_diameter * mesh.bbox.Diag());
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
        }
        
        if (verbose) std::cout << "[Writing file] ";
        filename.append(".");
        filename.append(format);
        if (format == "ply") {
            vcg::tri::io::ExporterPLY<TMesh>::Save(mesh, filename.c_str(), false);
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