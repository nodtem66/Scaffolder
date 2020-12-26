#include <string>
#include <sstream>
#include <algorithm>

#include "ProgressBar.hpp"
#include "utils.h"
#include "MeshOperation.h"
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

void print_menu() {
    std::cout << "(1) Fix self-intersect faces & border edges" << std::endl
        << "(2) Fix non-manifold vertices" << std::endl
        << "(3) Fix non-manifold edges" << std::endl
        << "(s) Save (export) and exit" << std::endl
        << "(e) Exit without save" << std::endl
        << "Menu: ";
}

int main(int argc, char* argv[])
{
    bool verbose = true;
    bool is_interaction_mode = false;
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
            ("i,interaction", "Interaction mode", cxxopts::value<bool>(), "BOOL")
            ("input", "Input file (STL,PLY,OFF,OBJ)", cxxopts::value<std::string>(), "FILE")
            ("output", "Output file (STL,PLY,OFF,OBJ)", cxxopts::value<std::string>(), "FILE")
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

        if (result.count("interaction")) {
            is_interaction_mode = true;
            verbose = true;
        }
    }
    catch (const cxxopts::OptionException & ex) {
        std::cout << "Error parsing options: " << ex.what() << std::endl;
        return 1;
    }

    std::ostream verbose_stream((verbose ? std::cout.rdbuf() : &util::null_buffer));
    
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

        clean_mesh(mesh, minimum_diameter, 0, verbose_stream);
        bool is_success = false;

        if (is_interaction_mode) {
            char menu;
            bool is_not_exit = true;
            do {
                report_mesh(verbose_stream, mesh);
                print_menu();
                std::cin >> menu;
                std::cout << std::endl;
                switch (menu) {
                case '1':
                    fix_self_intersect_mesh(mesh, minimum_diameter, 1, verbose_stream);
                    break;
                case '2':
                    fix_non_manifold_vertices(mesh, minimum_diameter, 1, verbose_stream);
                    break;
                case '3':
                    fix_non_manifold_vertices(mesh, minimum_diameter, 1, verbose_stream);
                    break;
                case 's':
                    is_success = true;
                case 'e':
                    is_not_exit = false;
                    break;
                }
                std::cout << std::endl;
            } while (is_not_exit);
            
        }
        else {
            // Automatic repairing
            is_success = fix_self_intersect_mesh(mesh, minimum_diameter, max_iteration, verbose_stream);

            bool is_manifold = is_mesh_manifold(mesh);
            if (!is_manifold) {
                is_success = fix_non_manifold_vertices(mesh, minimum_diameter, max_iteration, verbose_stream);
                is_success = fix_non_manifold_edges(mesh, minimum_diameter, max_iteration, verbose_stream);
            }
            if (verbose) report_mesh(verbose_stream, mesh);
        }

        

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
        else if (!is_success && !is_interaction_mode) {
            if (verbose) std::cout << "Sorry, fail to fix" << std::endl;
        }
    }
    return 0;
}