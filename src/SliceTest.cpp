#include <Eigen/Core>
#include <string>
#include <sstream>
#include <algorithm>
#include <ctime>

#include "cxxopts.hpp"
#include "ProgressBar.hpp"
#include "OptimalSlice.hpp"
#include "utils.h"

int main(int argc, char* argv[])
{
    // Define default parameters
    Eigen::IOFormat CleanFmt(4, Eigen::DontAlignCols, ", ", "\n", "[", "]");
    Eigen::IOFormat CSVFmt(-1, Eigen::DontAlignCols, ", ", ", ");
    // file parameters
    std::string input_file = "";
    size_t slice = 10;
    int direction = 2;
    bool is_export_convexhull = false;
    try {
        cxxopts::Options options("SliceTest", "Test");
        options.positional_help("[option args]").show_positional_help();
        options.add_options()
            ("h,help", "Print help")
            ("i,input", "Input file (STL)", cxxopts::value<std::string>(), "FILE")
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
        if (result.count("input")) input_file = result["input"].as<std::string>();
        else {
            std::cout << "Missing Input file" << std::endl;
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
    //slice::Polygon p = { {0,0}, {10,0}, {20,10}, {21, 20}, {10, 30}, {4, 31}, {-10, 22}, {-1, 20}, {8,18}, {2,6}, {-5,5} };
    //slice::Polygon p2 = { {-20,0}, {0,-20}, {20,0}, {0,20} };
    
    // Stage 1:
    {
        try {
            // Read FROM STL
            TMesh stl;
            int loadmark = 0;
            vcg::tri::io::ImporterSTL<TMesh>::Open(stl, input_file.c_str(), loadmark);
            std::cout << "STL Imported" << std::endl;
            vcg::tri::UpdateBounding<TMesh>::Box(stl);
            vcg::Box3d bbox = stl.bbox;
            vcg::Point3d dim = bbox.Dim();
            int next_direction = ((direction + 1) % 3);
            int prev_direction = ((direction + 2) % 3);

            vcg::tri::Clean<TMesh>::RemoveDuplicateVertex(stl);
            vcg::tri::Allocator<TMesh>::CompactEveryVector(stl);
            vcg::tri::UpdateTopology<TMesh>::FaceFace(stl);
            vcg::tri::Clean<TMesh>::RemoveDuplicateFace(stl);
            vcg::tri::UpdateTopology<TMesh>::FaceFace(stl);
            vcg::tri::Clean<TMesh>::RemoveZeroAreaFace(stl);
            vcg::tri::UpdateTopology<TMesh>::FaceFace(stl);
            std::cout << "Removed duplicate face and vertex" << std::endl;
            std::cout << "Update topology" << std::endl;
            vcg::tri::Clean<TMesh>::MergeCloseVertex(stl, SLICE_PRECISION*10);
            std::cout << "Merge close vertex" << std::endl;
            vcg::tri::UpdateTopology<TMesh>::FaceFace(stl);
            std::cout << "Update topology" << std::endl;
            
            int connectedComponentsNum = vcg::tri::Clean<TMesh>::CountConnectedComponents(stl);
            std::cout << "Mesh is composed by " << connectedComponentsNum << " connected component(s)" << std::endl;

            int edgeNum = 0, edgeBorderNum = 0, edgeNonManifoldNum = 0;
            vcg::tri::Clean<TMesh>::CountEdgeNum(stl, edgeNum, edgeBorderNum, edgeNonManifoldNum);
            int vertManifNum = vcg::tri::Clean<TMesh>::CountNonManifoldVertexFF(stl, true);
            if (edgeNonManifoldNum == 0 && vertManifNum == 0) {
                std::cout << "Mesh is two-manifold " << std::endl;
                int holeNum = vcg::tri::Clean<TMesh>::CountHoles(stl);
                std::cout << "Mesh has " << holeNum << " holes" << std::endl;
                int genus = vcg::tri::Clean<TMesh>::MeshGenus(stl.vn, edgeNum, stl.fn, holeNum, connectedComponentsNum);
                std::cout << "Genus is " << genus << std::endl;
            }

            slice::Slice s = slice::incremental_slicing(stl, slice, direction);
            slice::ContourSlice C = slice::contour_construct(s, direction);
            
            // create output directory
            if (!C.empty()) {
                size_t firstindex = input_file.find_last_of("/\\");
                firstindex = firstindex == string::npos ? 0 : firstindex + 1;
                size_t lastindex = input_file.find_last_of(".");
                std::string dir = input_file.substr(firstindex, lastindex - firstindex) + "_svg";
                std::cout << "Creating direction: " << dir << std::endl;
                make_dir(dir);

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
                slice::measure_feret_and_shape(C, minFeret, maxFeret, shape);
                std::copy(minFeret.begin(), minFeret.end(), std::ostream_iterator<double>(std::cout, " "));
                std::cout << std::endl;
                std::copy(maxFeret.begin(), maxFeret.end(), std::ostream_iterator<double>(std::cout, " "));
                std::cout << std::endl;
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