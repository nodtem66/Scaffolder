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
    try {
        cxxopts::Options options("SliceTest", "Test");
        options.positional_help("[option args]").show_positional_help();
        options.add_options()
            ("h,help", "Print help")
            ("i,input", "Input file (STL)", cxxopts::value<std::string>(), "FILE")
            ("s,slice", "Number of slice > 1 (Default: 10)", cxxopts::value<std::size_t>(), "POSITVE INTEGER MORE THAN 1")
            ("d,direction", "Direction 0=X, 1=Y, 2=Z default: 2)", cxxopts::value<int>(), "{0,1,2}");
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
    }
    catch (const cxxopts::OptionException & ex) {
        std::cout << "Error parsing options: " << ex.what() << std::endl;
        return 1;
    }

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

            vcg::tri::Clean<TMesh>::RemoveDuplicateVertex(stl);
            vcg::tri::Allocator<TMesh>::CompactEveryVector(stl);
            std::cout << "Removed duplicate face and vertex" << std::endl;
            vcg::tri::UpdateTopology<TMesh>::FaceFace(stl);
            std::cout << "Update topology" << std::endl;
            vcg::tri::Clean<TMesh>::MergeCloseVertex(stl, SLICE_PRECISION*100);
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
            slice::ContourSlice C = slice::contour_construct(s);
            
            // create output directory
            if (!C.empty()) {
                size_t firstindex = input_file.find_last_of("/\\");
                firstindex = firstindex == string::npos ? 0 : firstindex + 1;
                size_t lastindex = input_file.find_last_of(".");
                std::string dir = input_file.substr(firstindex, lastindex - firstindex) + "_svg";
                std::cout << "Creating direction: " << dir << std::endl;
                make_dir(dir);

                for (slice::ContourSlice::const_iterator cs = C.begin(); cs != C.end(); cs++) {
                    size_t index = (size_t)(cs - C.begin() + 1);
                    std::cout << "Contour slice " << index << " size: " << cs->size() << " contour(s)" << std::endl;
                    std::stringstream name(dir);
                    name << dir << "/" << index << ".svg";
                    slice::write_svg(name.str(), *cs, dim[0], dim[1], bbox.min[0], bbox.min[1]);
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