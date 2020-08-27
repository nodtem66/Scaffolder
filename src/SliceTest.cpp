#include <string>
#include <sstream>
#include <algorithm>
#include <ctime>

#include "cxxopts.hpp"
#include "ProgressBar.hpp"
#include "OptimalSlice.hpp"
#include "utils.h"
#include "MeshOperation.h"

ProgressBar import_progress(100, 40);
ProgressBar qsim_progress(100, 40);
int _pos = 0;
bool import_callback(int pos, const char* str) {
    if (pos % 10 == 0 && _pos != pos / 10) {
        _pos = pos / 10;
        import_progress.update(pos);
        import_progress.display();
    }
    if (pos >= 100) import_progress.done();
    return true;
}

int main(int argc, char* argv[])
{
    // file parameters
    std::string input_file = "";
    std::string format = "";
    size_t slice = 10;
    int direction = 2;
    bool verbose = true;
    bool is_export_convexhull = false;
    try {
        cxxopts::Options options("SliceTest", "Test");
        options.positional_help("[option args]").show_positional_help();
        options.add_options()
            ("h,help", "Print help")
            ("i,input", "Input file (STL,PLY,OFF,OBJ)", cxxopts::value<std::string>(), "FILE")
            ("s,slice", "Number of slice > 1 (Default: 10)", cxxopts::value<std::size_t>(), "POSITVE INTEGER MORE THAN 1")
            ("d,direction", "Direction 0=X, 1=Y, 2=Z default: 2)", cxxopts::value<int>(), "{0,1,2}")
            ("x,convexhull", "Export convexhull to SVG", cxxopts::value<bool>());
        options.parse_positional({ "input", "slice", "direction" });
        bool isEmptyOption = (argc == 1);
        cxxopts::ParseResult result = options.parse(argc, argv);
        if (isEmptyOption || result.count("help")) {
            std::cout << options.help() << std::endl;
            return 0;
        }
        // Requirment
        if (result.count("input")) {
            input_file = result["input"].as<std::string>();
            std::string ext = util::PathGetExtension(input_file);
            if (!ext.empty()) {
                format = ext.substr(1);
                util::to_lower(format);
            }
        }
        else {
            std::cout << "Missing Input file" << std::endl;
            return 1;
        }

        if (format != "ply" && format != "obj" && format != "stl" && format != "off") {
            std::cout << "Invalid input extension: " << format << std::endl;
            return 1;
        }

        if (result.count("slice")) slice = result["slice"].as<std::size_t>();
        if (slice < 2) {
            std::cout << "Invalid slice: " << slice << std::endl;
        }
        if (result.count("direction")) direction = result["direction"].as<int>();
        if (direction < 0 && direction > 2) {
            std::cout << "Invalid direction: " << direction << std::endl;
        }
        if (result.count("convexhull")) is_export_convexhull = result["convexhull"].as<bool>();
    }
    catch (const cxxopts::OptionException & ex) {
        std::cout << "Error parsing options: " << ex.what() << std::endl;
        return 1;
    }

    {
        try {
            // Read FROM STL
            TMesh mesh;
            {
                if (verbose) std::cout << "Import " << format << ": " << input_file << std::endl;
                int loadmark = 0;
                if (format == "stl") {
                    vcg::tri::io::ImporterSTL<TMesh>::Open(mesh, input_file.c_str(), loadmark, import_callback);
                }
                else if (format == "ply") {
                    vcg::tri::io::ImporterPLY<TMesh>::Open(mesh, input_file.c_str(), import_callback);
                }
                else if (format == "off") {
                    vcg::tri::io::ImporterOFF<TMesh>::Open(mesh, input_file.c_str(), import_callback);
                }
                else if (format == "obj") {
                    vcg::tri::io::ImporterOBJ<TMesh>::Open(mesh, input_file.c_str(), loadmark, import_callback);
                }
                if (verbose) std::cout << std::endl;
            }
            std::cout << "Mesh Imported" << std::endl;
            vcg::tri::UpdateBounding<TMesh>::Box(mesh);
            
            vcg::Box3d bbox = mesh.bbox;
            vcg::Point3d dim = bbox.Dim();
            int next_direction = ((direction + 1) % 3);
            int prev_direction = ((direction + 2) % 3);

            vcg::tri::Clean<TMesh>::RemoveDuplicateVertex(mesh);
            vcg::tri::Allocator<TMesh>::CompactEveryVector(mesh);
            vcg::tri::UpdateTopology<TMesh>::FaceFace(mesh);
            vcg::tri::Clean<TMesh>::RemoveDuplicateFace(mesh);
            vcg::tri::UpdateTopology<TMesh>::FaceFace(mesh);
            vcg::tri::Clean<TMesh>::RemoveZeroAreaFace(mesh);
            vcg::tri::UpdateTopology<TMesh>::FaceFace(mesh);
            std::cout << "Removed duplicate face and vertex" << std::endl;
            std::cout << "Update topology" << std::endl;
            vcg::tri::Clean<TMesh>::MergeCloseVertex(mesh, SLICE_PRECISION*10);
            std::cout << "Merge close vertex" << std::endl;
            vcg::tri::UpdateTopology<TMesh>::FaceFace(mesh);
            std::cout << "Update topology" << std::endl;
            
            int connectedComponentsNum = vcg::tri::Clean<TMesh>::CountConnectedComponents(mesh);
            std::cout << "Mesh is composed by " << connectedComponentsNum << " connected component(s)" << std::endl;

            int edgeNum = 0, edgeBorderNum = 0, edgeNonManifoldNum = 0;
            vcg::tri::Clean<TMesh>::CountEdgeNum(mesh, edgeNum, edgeBorderNum, edgeNonManifoldNum);
            int vertManifNum = vcg::tri::Clean<TMesh>::CountNonManifoldVertexFF(mesh, true);
            if (edgeNonManifoldNum == 0 && vertManifNum == 0) {
                std::cout << "Mesh is two-manifold " << std::endl;
                int holeNum = vcg::tri::Clean<TMesh>::CountHoles(mesh);
                std::cout << "Mesh has " << holeNum << " holes" << std::endl;
                int genus = vcg::tri::Clean<TMesh>::MeshGenus(mesh.vn, edgeNum, mesh.fn, holeNum, connectedComponentsNum);
                std::cout << "Genus is " << genus << std::endl;
            }

            slice::Slice s = slice::incremental_slicing(mesh, slice, direction);
            slice::ContourSlice C = slice::contour_construct(s, direction);
            
            // create output directory
            if (!C.empty()) {
                size_t firstindex = input_file.find_last_of("/\\");
                firstindex = firstindex == string::npos ? 0 : firstindex + 1;
                size_t lastindex = input_file.find_last_of(".");
                std::string dir = input_file.substr(firstindex, lastindex - firstindex);
                dir.append("_svg");
                std::cout << "Creating direction: " << dir << std::endl;
                util::make_dir(dir);

#ifndef USE_PARALLEL
                for (slice::ContourSlice::const_iterator cs = C.begin(); cs != C.end(); cs++) {
                    size_t index = (size_t)(cs - C.begin() + 1);
                    std::stringstream name(dir);
                    name << dir << "/" << index << ".svg";
                    slice::write_svg(name.str(), *cs, dim[next_direction], dim[prev_direction], bbox.min[next_direction], bbox.min[prev_direction], is_export_convexhull);
                }
#else
                tbb::parallel_for(
                    tbb::blocked_range<size_t>(0, C.size()),
                    [&](const tbb::blocked_range<size_t> &r) {
                        for (size_t index = r.begin(); index < r.end(); index++) {
                            std::stringstream name(dir);
                            name << dir << "/" << (index) << ".svg";
                            slice::write_svg(name.str(), C[index], dim[next_direction], dim[prev_direction], bbox.min[next_direction], bbox.min[prev_direction], is_export_convexhull);
                        }
                    }
                );
#endif
                std::vector<double> minFeret, maxFeret, shape[5];
                std::ofstream result;
                slice::measure_feret_and_shape(C, minFeret, maxFeret, shape);
                if (minFeret.size() > 0) {
                    std::sort(minFeret.begin(), minFeret.end());
                    std::cout << "Min Feret (Min, Q1, Q2, Q3, Max) = (" << minFeret.at(0) << ',' << minFeret.at(minFeret.size() * 0.25) << ',' << minFeret.at(minFeret.size() / 2) << ',' << minFeret.at(minFeret.size() * 0.75) << ',' << minFeret.at(minFeret.size() - 1) << ')' << std::endl;
                    result.open(dir + "/" + "min_feret.csv", std::ofstream::out);
                    for (std::vector<double>::iterator it = minFeret.begin(); it != minFeret.end(); ++it) {
                        result << *it << ",";
                    }
                    result.close();
                }
                
                if (maxFeret.size() > 0) {
                    std::sort(maxFeret.begin(), maxFeret.end());
                    std::cout << "Max Feret (Min, Q1, Q2, Q3, Max) = (" << maxFeret.at(0) << ',' << maxFeret.at(maxFeret.size() * 0.25) << ',' << maxFeret.at(maxFeret.size() / 2) << ',' << maxFeret.at(maxFeret.size() * 0.75) << ',' << maxFeret.at(maxFeret.size() - 1) << ')' << std::endl;
                    result.open(dir + "/" + "max_feret.csv", std::ofstream::out);
                    for (std::vector<double>::iterator it = maxFeret.begin(); it != maxFeret.end(); ++it) {
                        result << *it << ",";
                    }
                    result.close();
                }
                
            }
            
        }
        catch (const std::exception & ex) {
            std::cout << "Exception: " << ex.what() << std::endl;
            std::cout << strerror(errno) << endl;
            return errno;
        }
    }
    return 0;
}