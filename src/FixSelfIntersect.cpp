#include <string>
#include <sstream>
#include <algorithm>

#include "ProgressBar.hpp"
#include "utils.h"
#include "Mesh.h"
#include "cxxopts.hpp"

ProgressBar import_progress(100, 40);

bool import_callback(int pos, const char* str) {
    if (pos % 10 == 0) {
        import_progress.update(pos);
        import_progress.display();
    }
    if (pos >= 100) import_progress.done();
    return true;
}

int main(int argc, char* argv[])
{
    bool verbose = true;
    double minimum_diameter = 0.25;
    uint16_t max_iteration = 5;

    std::string input_file = "";
    std::string format = "";
    std::string output_file = "";
    std::string output_format = "";

    try {
        cxxopts::Options options("Fix Self Intersect", "Automatic fix self-intersect");
        options.positional_help("INPUT OUTPUT MAX_ITERATION").show_positional_help();
        options.add_options()
            ("h,help", "Print help")
            ("i,input", "Input file (STL,PLY,OFF,OBJ)", cxxopts::value<std::string>(), "FILE")
            ("o,output", "Output file (STL,PLY,OFF,OBJ)", cxxopts::value<std::string>(), "FILE")
            ("n,max_iteration", "Max iteration [default: 5]", cxxopts::value<uint16_t>(), "UINT16");
        options.parse_positional({ "input", "output", "max_iteration" });
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

        // Optional parameter: output_file default to input_file
        if (result.count("output")) {
            output_file = result["output"].as<std::string>();
            std::string ext = util::PathGetExtension(output_file);
            if (!ext.empty()) {
                output_format = ext.substr(1);
                util::to_lower(output_format);
            }
        }
        else {
            output_file = input_file;
            output_format = format;
        }

        if (format != "ply" && format != "obj" && format != "stl" && format != "off") {
            std::cout << "Invalid input extension: " << format << std::endl;
            return 1;
        }

        if (output_format != "ply" && output_format != "obj" && output_format != "stl" && output_format != "off") {
            std::cout << "Invalid output extension: " << output_format << std::endl;
            return 1;
        }

        if (result.count("max_iteration")) max_iteration = result["max_iteration"].as<uint16_t>();
    }
    catch (const cxxopts::OptionException & ex) {
        std::cout << "Error parsing options: " << ex.what() << std::endl;
        return 1;
    }
    
    {
        TMesh mesh;
        {
            if (verbose) std::cout << "Import " << format << std::endl;
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
        vcg::tri::UpdateBounding<TMesh>::Box(mesh);

        bool is_success = fix_self_intersect_mesh(mesh, minimum_diameter, max_iteration, verbose);

        bool is_manifold = is_mesh_manifold(mesh);
        if (!is_manifold) {
            is_success = fix_non_manifold(mesh, minimum_diameter, max_iteration, verbose);
        }
        if(verbose) report_mesh(mesh);

        if (is_success) {
            if (verbose) std::cout << "Export " << output_file;
            if (output_format == "stl") {
                vcg::tri::io::ExporterSTL<TMesh>::Save(mesh, output_file.c_str(), true);
            }
            else if (output_format == "ply") {
                vcg::tri::io::ExporterPLY<TMesh>::Save(mesh, output_file.c_str(), false);
            }
            else if (output_format == "off") {
                vcg::tri::io::ExporterOFF<TMesh>::Save(mesh, output_file.c_str());
            }
            else if (output_format == "obj") {
                vcg::tri::io::ExporterOBJ<TMesh>::Save(mesh, output_file.c_str(), 0);
            }
            if (verbose) std::cout << " [OK]" << std::endl;
        }
        else {
            if (verbose) std::cout << "Sorry, fail to fix" << std::endl;
        }
    }
    return 0;
}