#include "Scaffolder_2.h"
#include "QuadricSimplification.h"

ProgressBar qsim_progress(100, 40);
inline bool qsim_callback(int pos, const char* str) {
    if (pos >= 0 && pos <= 100) {
        qsim_progress.update(pos);
        qsim_progress.display();
    }
    if (pos >= 100)
        qsim_progress.done();
    return true;
}

int main(int argc, char* argv[])
{
    // Define default parameters
    bool verbose = true;
    bool dirty = false;
    bool is_analysis_microstructure = false;
    bool is_export_microstructure = false;
    bool is_mean_curvature = false;
    bool is_export_feret = false;
    bool is_build_inverse = false;
    bool is_fix_self_intersect = false;
    bool is_export_jpeg = false;
    bool no_output = true;

    uint8_t export_axis = 'X';
    uint16_t grid_offset = 3;
    uint16_t shell = 0;
    uint16_t smooth_step = 5;
    uint16_t k_slice = 100;
    uint16_t k_polygon = 4;
    size_t minimum_grid_size = 100;

    double isolevel = 0.0;
    double qsim_percent = 0;
    double coff = pi;
    double minimum_diameter = 0.25;

    // format Eigen object to readable text
    Eigen::IOFormat CleanFmt(4, Eigen::DontAlignCols, ", ", "\n", "[", "]");
    Eigen::IOFormat CSVFmt(-1, Eigen::DontAlignCols, ", ", ", ");

    // file parameters
    std::string filename = "", dir = "";
    std::string surface_name = "schwarzp";
    std::string input_file = "";
    std::string output_format = "";
    std::uint8_t log_format = SCAFFOLDER_FORMAT_DEFAULT;

    // TRY-CATCH block
    // Parse the program options
    try {
        cxxopts::Options options("Scaffolder", "Scaffolder - generate 3D scaffold from STL file based on implicit surface");
        options.positional_help("INPUT OUTPUT PARAMETERS").show_positional_help();
        options.add_options()
            ("h,help", "Print help")
            ("i,input", "Input file (STL/PLY/OFF/OBJ/VMI)", cxxopts::value<std::string>(), "INPUT")
            ("o,output", "Output filename with extension stl,ply,obj,off,ctm", cxxopts::value<std::string>(), "OUTPUT")
            ("params", "Combined parameters list: surface[,coff,isolevel,grid_size,k_slice,k_polygon]", cxxopts::value<std::vector<std::string>>(), "PARAMETERS")
            ("q,quiet", "Disable verbose output [default: false]")
            ("c,coff", "Angular frequency (pore size adjustment) default:PI", cxxopts::value<double>(), "DOUBLE")
            ("t,isolevel", "isolevel (porosity adjustment) [default: 0]", cxxopts::value<double>(), "DOUBLE")
            ("n,surface", "implicit surface: rectlinear, schwarzp, schwarzd, gyroid, double-p, double-d, double-gyroiod, lidinoid, schoen_iwp, neovius, bcc, tubular_g_ab, tubular_g_c [default: schwarzp]", cxxopts::value<std::string>(), "NAME")
            ("g,grid_size", "Grid size [default: 100]", cxxopts::value<size_t>(), "INT (0..60000)")
            ("s,shell", "Outer thickness (layers) [default:0]", cxxopts::value<uint16_t>(), "INT (0..60000)")
            ("grid_offset", "[default:3]", cxxopts::value<uint16_t>(), "INT (0..60000)")
            ("m,microstructure", "Analysis microstructure with Slice contour technique ( [default: false]")
            ("export_microstructure", "Analysis microstructure and export the 2D contours (for debugging) [default: false]")
            ("export_jpeg", "Export 2D JPEG [Default: Z,100]", cxxopts::value<std::vector<std::string>>(), "[X|Y|Z],INT")
            ("k_slice", "K_slice: the number of slicing layers in each direction (used in microstructure analysis) [default: 100]", cxxopts::value<uint16_t>(), "INT (0..60000)")
            ("k_polygon", "K_polygon: the number of closest outer contour (used in microstructure analysis) [default: 4]", cxxopts::value<uint16_t>(), "INT (>0)")
            ("z,size_optimize", "Experimental Quadric simplification [default: 0]", cxxopts::value<double>(), "DOUBLE (0..1)")
            ("smooth_step", "Smooth with laplacian (default: 5)", cxxopts::value<uint16_t>(), "INT (0..60000)")
            ("dirty", "Disable autoclean [default false]")
            ("minimum_diameter", "used for removing small orphaned (between 0-1) [default: 0.25]", cxxopts::value<double>(), "DOUBLE (0..1)")
            ("format", "Format of logging output [default: default]", cxxopts::value<std::string>(), "FORMAT (default, csv)")
            ("output_inverse", "additional output inverse scaffold [default: false]")
            ("fix_self_intersect", "Experimental fix self-intersect faces [default: false]")
            ("mean_curvature", "Experimental calculate mean curvature");
        options.parse_positional({ "input", "output", "params" });
        bool isEmptyOption = (argc == 1);
        cxxopts::ParseResult result = options.parse(argc, argv);
        if (isEmptyOption || result.count("help")) {
            std::cout << options.help() << std::endl
                << "Example: " << std::endl << std::endl
                << "  " << util::PathGetBasename(argv[0]) << " input.stl output.stl bcc,3.14159,0,100" << std::endl
                << "    " << "Generated BCC scaffold with w=3.14159 (PI), t=0, and grid size=100" << std::endl << std::endl
                << "  " << util::PathGetBasename(argv[0]) << " input.stl output.stl custom.lua,3.14159,0,100,100,4 -m" << std::endl
                << "    " << "Generated and evaluated scaffold with custom.lua, w=3.14159 (PI), t=0, " << std::endl
                << "    " << "grid size=100, k_slice=100, k_polygon=4" << std::endl << std::endl
                << "  " << util::PathGetBasename(argv[0]) << " input.stl output.stl bcc,3.14159,0 -m -q --format csv" << std::endl
                << "    " << "Generated and evaluated BCC scaffold (w=3.14159, t=0) and reported in CSV" << std::endl << std::endl << std::endl
                << "Lua file: " << std::endl << std::endl
                << "  Define the \"surface\" function that return the implicit function" << std::endl
                << "  -----------------------------------------------------------------" << std::endl
                << "  function surface (x, y, z)" << std::endl
                << "    return sin(x) * cos(y) + sin(y) * cos(z) + sin(z) * cos(x) - params.isolevel" << std::endl
                << "  end" << std::endl
                << "  -----------------------------------------------------------------" << std::endl << std::endl
                << "Special symbols can be used in lua file: " << std::endl << std::endl
                << "  params = { coff, isolevel, k_splice, k_polygon }" << std::endl
                << "  bbox = { min, max, length, grid_density }" << std::endl
                << "  winding(x,y,z): function returning the winding number of position x,y,z" << std::endl
                << "  signed_distance(x,y,z): function returning signed distance of position x,y,z" << std::endl
                << "  and all functions from math module" << std::endl << std::endl;
            return 0;
        }
        // Program requires at least one argument for specifying INPUT file
        // Show error if there's no argument
        if (result.count("input")) input_file = result["input"].as<std::string>();
        else {
            std::cout << "Missing Input file" << std::endl;
            return 1;
        }
        // Parse the optional parameters
        if (result.count("quiet")) verbose = !result["quiet"].as<bool>();
        if (result.count("format")) {
            std::string fmt = result["format"].as<std::string>();
            util::to_lower(fmt);
            if (fmt == "csv") {
                log_format = SCAFFOLDER_FORMAT_CSV;
            }
        }
        if (result.count("dirty")) dirty = result["dirty"].as<bool>();
        if (result.count("microstructure")) {
            is_analysis_microstructure = result["microstructure"].as<bool>();
            no_output = true;
        }
        if (result.count("export_microstructure")) {
            is_analysis_microstructure = result["m2"].as<bool>();
            is_export_microstructure = result["m2"].as<bool>();
            no_output = true;
        }
        if (result.count("output_inverse")) {
            is_build_inverse = result["output_inverse"].as<bool>();
        }
        if (result.count("output")) {
            filename = result["output"].as<std::string>();
            util::to_lower(filename);
            std::string ext = util::PathGetExtension(filename);
            if (!ext.empty()) {
                filename = filename.substr(0, filename.size() - ext.size());
                output_format = ext.substr(1);
                if (output_format != "ply" && output_format != "obj" && output_format != "stl" && output_format != "off" && output_format != "ctm") {
                    std::cout << "Invalid format: " << output_format << std::endl;
                    return 1;
                }
                no_output = false;
            }
        }
        if (result.count("params")) {
            std::vector<std::string> params = result["params"].as<std::vector<std::string>>();
            size_t n = params.size();
            if (n >= 1) surface_name = params[0];
            if (n >= 2) coff = std::stod(params[1]);
            if (n >= 3) isolevel = std::stod(params[2]);
            if (n >= 4) minimum_grid_size = std::stoi(params[3]);
            if (n >= 5) k_slice = std::stoi(params[4]);
            if (n >= 6) k_polygon = std::stoi(params[5]);
        }
        if (result.count("export_jpeg")) {
            std::vector<std::string> params = result["export_jpeg"].as<std::vector<std::string>>();
            size_t n = params.size();
            is_export_jpeg = true;
            if (n >= 1) export_axis = toupper(params[0][0]);
            if (n >= 2) k_slice = std::stoi(params[1]);
        }
        if (result.count("isolevel")) isolevel = result["isolevel"].as<double>();
        if (result.count("grid_size")) minimum_grid_size = result["grid_size"].as<size_t>();
        if (result.count("grid_offset")) grid_offset = result["grid_offset"].as<uint16_t>();
        if (result.count("coff")) coff = result["coff"].as<double>();
        if (result.count("minimum_diameter")) minimum_diameter = result["minimum_diameter"].as<double>();
        if (result.count("surface")) surface_name = result["surface"].as<std::string>();
        if (result.count("shell")) shell = result["shell"].as<uint16_t>();
        if (result.count("smooth_step")) smooth_step = result["smooth_step"].as<uint16_t>();
        if (result.count("fix_self_intersect")) is_fix_self_intersect = result["fix_self_intersect"].as<bool>();
        if (result.count("size_optimize")) qsim_percent = result["size_optimize"].as<double>();
        if (result.count("k_slice")) k_slice = result["k_slice"].as<uint16_t>();
        if (result.count("k_polygon")) k_polygon = result["k_polygon"].as<uint16_t>();
        if (result.count("mean_curvature")) is_mean_curvature = result["mean_curvature"].as<bool>();

        util::to_lower(surface_name);
        util::to_lower(output_format);

        if (surface_name == "rectlinear") {
            isolevel = 0;
        }

        if (is_export_microstructure || is_export_jpeg) {
            // Initialize filename and directory for logging
            // Get filename without extension from input_file.
            // Ex., input.stl -> input
            // then the folder `input_slice` will be created
            dir = util::PathGetBasename(input_file) + "_slice";
            util::make_dir(dir);
        }
    }
    catch (const cxxopts::OptionException& ex) {
        std::cout << "Error parsing options: " << ex.what() << std::endl;
        return 1;
    }

    // Initialize output method via log file or standard output
    std::ofstream log_file;
    std::streambuf* log_buffer;
    time_t run_timestamp = std::time(nullptr);
    if (verbose) {
        log_buffer = std::cout.rdbuf();
    }
    else {
        std::stringstream _name;
        _name << filename << '_' << surface_name << run_timestamp << "." << ((log_format == SCAFFOLDER_FORMAT_CSV) ? "csv" : "txt");
        log_file.open(_name.str(), std::ofstream::out);
        log_buffer = log_file.rdbuf();
    }
    std::ostream log(log_buffer);
    util::NullBuffer null_buffer;
    std::ostream verbose_stream((verbose ? log_buffer : &null_buffer));

    LOG << "[Scaffolder " << VERSION << "]" << std::endl
        << "-- Input file: " << input_file << std::endl
        << "-- Output file: " << filename << '.' << output_format << std::endl
        << "-- Surface (-n): " << surface_name << std::endl
        << "-- Coff (-c): " << coff << std::endl
        << "-- Isolevel (-t): " << isolevel << std::endl
        << "-- Grid size (-g): " << minimum_grid_size << std::endl
        << "--   Grid offset: " << grid_offset << std::endl
        << "--   Shell: " << shell << std::endl
        << "-- Autoclean: " << (dirty ? "False" : "True") << std::endl
        << "--   Minimum diameter: " << 100 * minimum_diameter << "%" << std::endl
        << "--   Smooth step: " << smooth_step << std::endl
        << "--   Fix self-intersect: " << (is_fix_self_intersect ? "True" : "False") << std::endl
        << "--   Quadric Simplification (-z): " << qsim_percent << std::endl
        << "-- Analysis microstructure (-m): " << (is_analysis_microstructure ? "True" : "False") << std::endl
        << "--   Slice grid (k_slice): " << k_slice << std::endl
        << "--   Nearest outer contours (k_polygon): " << k_polygon << std::endl
        << "--   Export microstructure: " << (is_export_microstructure ? "True" : "False") << std::endl
        << "--   Mean curvature: " << (is_mean_curvature ? "True" : "False") << std::endl
        << "-- Export JPEG: " << (is_export_jpeg ? "Yes" : "No") << std::endl
        << "--   Axis: " << export_axis << std::endl
        << "-- Build: " << (no_output ? "No" : "Yes") << (is_build_inverse ? " (Inverse)" : "") << std::endl;
    else
    CSV << "Surface,coff,shell,thickness,grid_size,grid_offset,smooth_step,input_file,";

    if (is_analysis_microstructure)
        CSV << "avg_min_feret,avg_max_feret,min_min_feret,q1_min_feret,q2_min_feret,q3_min_feret,max_min_feret,"
            << "min_max_feret,q1_max_feret,q2_max_feret,q3_max_feret,max_max_feret,"
            << "min_square,q1_square,q2_square,q3_square,max_square,"
            << "min_circle,q1_circle,q2_circle,q3_circle,max_circle,"
            << "min_triangle,q1_triangle,q2_triangle,q3_triangle,max_triangle,"
            << "min_ellipse,q1_ellipse,q2_ellipse,q3_ellipse,max_ellipse,"
            << "min_elongation,q1_elongation,q2_elongation,q3_elongation,max_elongation,";
    if (is_mean_curvature)
        CSV << "min_mean_curvature,q1_mean_curvature,q2_mean_curvature,q3_mean_curvature,max_mean_curvature,";

    CSV << "volumn,surface_area,porosity,surface_area_ratio,vertices,faces" << std::endl
        << surface_name << ',' << coff << ',' << shell << ',' << isolevel << ',' << minimum_grid_size << ',' << grid_offset << ',' << smooth_step << ','
        << input_file << ',';

    // Next, we divided into two stages separated by scope because of the need of memory cleaning of unused variables
    // Stage 1
    TMesh mesh, inverse_mesh;
    // volumn1 and area1 is volumn and surface area of original mesh
    // volumn2 and area2 is volumn and surface area of generated mesh
    // we initialize these variables to a small number (epsilon) because of prevention of divided-by-zero exception
    double volume1 = eps, volume2 = eps;
    double area1 = eps, area2 = eps;
    Eigen::RowVector3d V1min1;
    Eigen::RowVector3i grid_size;
    double grid_delta;
    bool is_manifold;
    try {
        // Read vertices and faces from STL
        Eigen::MatrixXd Fxyz, IFxyz;
        Eigen::MatrixXd V1;
        Eigen::MatrixXi F1;
        Eigen::RowVector3d V1min, V1max, L, delta;
        {
            TMesh stl;
            int loadmark = 0;
            int err = vcg::tri::io::Importer<TMesh>::Open(stl, input_file.c_str());
            if (err)
            {
                std::cout << "Unable to open mesh " << input_file << " : " << vcg::tri::io::Importer<TMesh>::ErrorMsg(err) << std::endl;
                exit(-1);
            }
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

            // Calculate the grid size parameters
            V1min = V1.colwise().minCoeff();
            V1max = V1.colwise().maxCoeff();
            L = V1max - V1min;
            delta = L / minimum_grid_size;
            grid_delta = delta.minCoeff();
            // Create border offset from the original box
            V1min1 = V1min - grid_offset * grid_delta * Eigen::RowVector3d::Ones();
            grid_size = (L / grid_delta).cast<int>() + 2 * grid_offset * Eigen::RowVector3i::Ones();
            LOG << "-- Bounding Box: "
                << V1min.format(CleanFmt) << ' ' << V1max.format(CleanFmt) << std::endl
                << "-- Length: " << L.format(CleanFmt) << std::endl
                << "-- Grid delta: " << grid_delta << std::endl;
        }

        // Generate Grid index (x,y,z)
        LOG << "[Generating grid] ";
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
        LOG << "OK" << std::endl
            << "-- Grid size: " << grid_size.prod() << " " << grid_size.format(CleanFmt) << std::endl
            << "[Calculating Winding number] ";
            
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
        LOG << "OK" << std::endl
            << "-- Sign Distance: [ " << Signed_Distance.minCoeff() << ", " << Signed_Distance.maxCoeff() << "]  "
            << "Wind: [ " << W.minCoeff() << ", " << W.maxCoeff() << "]" << std::endl;
            
        LOG << "[Generating isosurface Fxyz] ";
        sol::state lua;
        Function* surface;
        sol::function f;

        // Initialize surface
        bool isLuaFunction = false;
        if (util::PathGetExtension(surface_name) == ".lua") {
            lua.open_libraries(sol::lib::base, sol::lib::math, sol::lib::string);
            set_shorten_function(lua);
            sol::load_result lua_file = lua.load_file(surface_name);
            if (!lua_file.valid()) {
                sol::error err = lua_file;
                std::cerr << "[Lua Error] " << err.what() << std::endl;
                exit(-1);
            }
            else {
                lua["params"] = lua.create_table_with(
                    "coff", coff,
                    "isolevel", isolevel,
                    "k_slice", k_slice,
                    "k_polygon", k_polygon
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
            surface = isosurface(surface_name, isolevel);
        }

        Implicit_function fn(surface, coff);
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
                    if (w < 0.5) {
                        Fxyz(index) = -1.0;
                        IFxyz(index) = -1.0;
                    }
                    else {
                        if (shell > 0 && s <= shell) {
                            Fxyz(index) = s;
                            IFxyz(index) = s;
                        }
                        else {
                            Fxyz(index) = fv;
                            IFxyz(index) = -fv;
                        }
                    }
                }
            }
        }
        LOG << "OK" << std::endl;

        // Export JPEG
        if (is_export_jpeg) {
            LOG << "[Export JPEG] ";
            const uint8_t Z = export_axis - 'X', X = (Z + 1) % 3, Y = (X + 1) % 3;
            const double scale = grid_size(Z) / k_slice;
            const uint16_t width = grid_size(X) - 2 * grid_offset;
            const uint16_t height = grid_size(Y) - 2 * grid_offset;
            for (uint16_t z_index = 0; z_index < k_slice; z_index++) {

                auto z = (uint16_t)floor(z_index * scale);

                std::stringstream name(dir);

                name << dir << "/" << export_axis << '_' << std::setfill('0') << std::setw(4) << (z) << ".jpg";
                std::ofstream jpg(name.str(), std::ios_base::out | std::ios_base::binary);

                auto image = new unsigned char[size_t(width) * height];

                for (auto y = 0; y < height; y++) {
                    for (auto x = 0; x < width; x++) {
                        const size_t index = (size_t(x)+grid_offset) + size_t(grid_size(X)) * ((y+grid_offset) + grid_size(Y) * (z+grid_offset));
                        auto offset = x + width * y;
                        image[offset] = (Fxyz(index) >= 0 ? 0 : 255);
                    }
                }
                TooJpeg::writeJpeg(jpg, image, width, height, false, 80, true);
                jpg.close();
                delete[] image;
            }
            std::ofstream infofile(dir + "/" + "info.txt", std::ios_base::out);
            infofile << "{" << std::endl
                << "\t\"unit per pixel\": " << grid_delta << std::endl
                << "\t\"width\": " << width << std::endl
                << "\t\"height\": " << height << std::endl
                << "}";
            infofile.close();
            LOG << "[OK]" << std::endl;
        } // End Export JPEG

        // Create scaffolder and inverse mesh by dual marching cube
        marching_cube(mesh, Fxyz, grid_size, V1min1, grid_delta, verbose, false);
        LOG << "-- Info: " << mesh.VN() << " vertices " << mesh.FN() << " faces" << std::endl;
        if (is_build_inverse) marching_cube(inverse_mesh, IFxyz, grid_size, V1min1, grid_delta, false, false);

        if (!dirty) {
            clean_mesh(mesh, minimum_diameter, smooth_step, verbose_stream);
            if (is_fix_self_intersect) fix_self_intersect_mesh(mesh, minimum_diameter, 5, verbose_stream);
            if (is_build_inverse) clean_mesh(inverse_mesh, minimum_diameter, 0, verbose_stream);
        }

        if (qsim_percent > 0) {
            LOG << "[Quadric Simplificaion: " << qsim_percent << "]" << std::endl;
            qsim_progress.reset();
            vcg::tri::mesh_quad_simplification(mesh, qsim_percent, qsim_callback);
            qsim_progress.done();
            LOG << "-- Info: " << mesh.VN() << " vertices " << mesh.FN() << " faces" << std::endl;
        }

        is_manifold = is_mesh_manifold(mesh);
        LOG << "-- is_manifold: " << is_manifold << std::endl;
        if (!is_manifold && !no_output) {
            is_manifold = fix_non_manifold_vertices(mesh, minimum_diameter, 5, verbose_stream);
            is_manifold = is_manifold && fix_non_manifold_edges(mesh, minimum_diameter, 5, verbose_stream);
            LOG << "-- is_manifold: " << is_manifold << std::endl;
        }
        if (!is_manifold && !verbose) {
            LOG << "[Warning] Mesh is not two manifold" << std::endl;
        }

        if (log_format == SCAFFOLDER_FORMAT_DEFAULT) report_mesh(log, mesh);
    }
    catch (const std::exception& ex) {
        std::cout << "Exception: " << ex.what() << std::endl;
        exit(-1);
    }

    try {
        if (is_analysis_microstructure && is_manifold) {
            // Evaluating pore size by create 2D slice from 3D mesh
            
            std::stringstream filename;
            std::vector<double> minFeret;
            std::vector<double> maxFeret;
            std::vector<double> podczeckShapes[5];
            
            // init progress bar
            ProgressBar progress(grid_size.sum(), PROGRESS_BAR_COLUMN);
                
            // Start measure pore size by Slice contour technique
            vcg::tri::UpdateBounding<TMesh>::Box(mesh);
            vcg::Box3f bbox = mesh.bbox;
            vcg::Point3f dim = bbox.Dim();
                    
            char axis = 'x';
            progress.display();
            // Loop for x, y, z axis
            // main_axis: 0 = x, 1 = y, 2 = z
            for (uint8_t main_axis = 0; main_axis < 3; main_axis++) {
                // define next and prev axis by the cycled direction: X -> Y, Y -> Z, Z -> X
                int next_axis = ((main_axis + 1) % 3);
                int prev_axis = ((main_axis + 2) % 3);
                // define progress bar step based on gize size
                unsigned int progress_major_step = grid_size[main_axis] / 4;
                // Rodrigo's incremental slicing
                optimal_slice::Slice s = optimal_slice::incremental_slicing(mesh, k_slice, main_axis);
                // update progress bar
                progress += progress_major_step;
                progress.display();
                // Rodrigo's contour construct
                optimal_slice::ContourSlice C = optimal_slice::contour_construct(s, main_axis);
                // update progress bar
                progress += progress_major_step;
                progress.display();
                // Measure pore size from C, and store the output in minFeret, maxFeret, and podczeckShapes
                optimal_slice::measure_feret_and_shape(C, k_polygon, minFeret, maxFeret, podczeckShapes);
                // update progress bar
                progress += progress_major_step;
                progress.display();
                // If export flag was defined, convert the 2D contours into SVG
                if (is_export_microstructure) {
                    filename.str(std::string());
                    unsigned int minor_step = C.size() / progress_major_step, step = minor_step;
                    // Loop foreach 2D contours
                    for (optimal_slice::ContourSlice::const_iterator cs = C.begin(); cs != C.end(); cs++) {
                        // Get array index from iterator
                        size_t index = (size_t)(cs - C.begin() + 1);
                        std::stringstream name(dir);
                        // Set the SVG filename to <axis>_<index>.svg
                        // Ex. x_0.svg, x_1.svg, ..., z_0.svg 
                        name << dir << "/" << axis << '_' << index << ".svg";
                        optimal_slice::write_svg(name.str(), *cs, dim[next_axis], dim[prev_axis], bbox.min[next_axis], bbox.min[prev_axis]);
                        step--;
                        if (step == 0) {
                            step = minor_step;
                            ++progress;
                            progress.display();
                        }
                    }
                }
                else {
                    progress += progress_major_step;
                    progress.display();
                }
                axis++;
            }
            // End Slice contour
            progress.done();

            // If we can measure pore sizes, then format the result
            if (minFeret.size() > 0 && maxFeret.size() > 0) {
                std::sort(minFeret.begin(), minFeret.end());
                std::sort(maxFeret.begin(), maxFeret.end());
                for (int i = 0; i < 5; i++) {
                    if (podczeckShapes[i].empty()) {
                        podczeckShapes[i].push_back(0);
                    }
                    else {
                        std::sort(podczeckShapes[i].begin(), podczeckShapes[i].end());
                    }
                }
                LOG << "[Microstructure] " << std::endl
                    << "-- Avg Min Feret: " << std::accumulate(minFeret.begin(), minFeret.end(), 0.0) / minFeret.size() << std::endl
                    << "-- Avg Max Feret: " << std::accumulate(maxFeret.begin(), maxFeret.end(), 0.0) / maxFeret.size() << std::endl
                    << "-- Min Feret: [" << minFeret.at(0) << " " << minFeret.at(minFeret.size() * 0.25) << " " << minFeret.at(minFeret.size() / 2) << " " << minFeret.at(minFeret.size() * 0.75) << " " << minFeret.at(minFeret.size() - 1) << "]" << std::endl
                    << "-- Max Feret: [" << maxFeret.at(0) << " " << maxFeret.at(maxFeret.size() * 0.25) << " " << maxFeret.at(maxFeret.size() / 2) << " " << maxFeret.at(maxFeret.size() * 0.75) << " " << maxFeret.at(maxFeret.size() - 1) << "]" << std::endl
                    << "-- Square Similarity: [" << podczeckShapes[0].at(0) << " " << podczeckShapes[0].at(podczeckShapes[0].size() * 0.25) << " " << podczeckShapes[0].at(podczeckShapes[0].size() / 2) << " " << podczeckShapes[0].at(podczeckShapes[0].size() * 0.75) << " " << podczeckShapes[0].at(podczeckShapes[0].size() - 1) << "]" << std::endl
                    << "-- Circle Similarity: [" << podczeckShapes[1].at(0) << " " << podczeckShapes[1].at(podczeckShapes[1].size() * 0.25) << " " << podczeckShapes[1].at(podczeckShapes[1].size() / 2) << " " << podczeckShapes[1].at(podczeckShapes[1].size() * 0.75) << " " << podczeckShapes[1].at(podczeckShapes[1].size() - 1) << "]" << std::endl
                    << "-- Triangle Similarity: [" << podczeckShapes[2].at(0) << " " << podczeckShapes[2].at(podczeckShapes[2].size() * 0.25) << " " << podczeckShapes[2].at(podczeckShapes[2].size() / 2) << " " << podczeckShapes[2].at(podczeckShapes[2].size() * 0.75) << " " << podczeckShapes[2].at(podczeckShapes[2].size() - 1) << "]" << std::endl
                    << "-- Ellipse Similarity: [" << podczeckShapes[3].at(0) << " " << podczeckShapes[3].at(podczeckShapes[3].size() * 0.25) << " " << podczeckShapes[3].at(podczeckShapes[3].size() / 2) << " " << podczeckShapes[3].at(podczeckShapes[3].size() * 0.75) << " " << podczeckShapes[3].at(podczeckShapes[3].size() - 1) << "]" << std::endl
                    << "-- Elongation Similarity: [" << podczeckShapes[4].at(0) << " " << podczeckShapes[4].at(podczeckShapes[4].size() * 0.25) << " " << podczeckShapes[4].at(podczeckShapes[4].size() / 2) << " " << podczeckShapes[4].at(podczeckShapes[4].size() * 0.75) << " " << podczeckShapes[4].at(podczeckShapes[4].size() - 1) << "]" << std::endl;
                else 
                CSV << std::accumulate(minFeret.begin(), minFeret.end(), 0.0) / minFeret.size() << ','
                    << std::accumulate(maxFeret.begin(), maxFeret.end(), 0.0) / maxFeret.size() << ','
                    << minFeret.at(0) << ',' << minFeret.at(minFeret.size() * 0.25) << ',' << minFeret.at(minFeret.size() / 2) << ',' << minFeret.at(minFeret.size() * 0.75) << ',' << minFeret.at(minFeret.size() - 1) << ','
                    << maxFeret.at(0) << ',' << maxFeret.at(maxFeret.size() * 0.25) << ',' << maxFeret.at(maxFeret.size() / 2) << ',' << maxFeret.at(maxFeret.size() * 0.75) << ',' << maxFeret.at(maxFeret.size() - 1) << ','
                    << podczeckShapes[0].at(0) << ',' << podczeckShapes[0].at(podczeckShapes[0].size() * 0.25) << ',' << podczeckShapes[0].at(podczeckShapes[0].size() / 2) << ',' << podczeckShapes[0].at(podczeckShapes[0].size() * 0.75) << ',' << podczeckShapes[0].at(podczeckShapes[0].size() - 1) << ','
                    << podczeckShapes[1].at(0) << ',' << podczeckShapes[1].at(podczeckShapes[1].size() * 0.25) << ',' << podczeckShapes[1].at(podczeckShapes[1].size() / 2) << ',' << podczeckShapes[1].at(podczeckShapes[1].size() * 0.75) << ',' << podczeckShapes[1].at(podczeckShapes[1].size() - 1) << ','
                    << podczeckShapes[2].at(0) << ',' << podczeckShapes[2].at(podczeckShapes[2].size() * 0.25) << ',' << podczeckShapes[2].at(podczeckShapes[2].size() / 2) << ',' << podczeckShapes[2].at(podczeckShapes[2].size() * 0.75) << ',' << podczeckShapes[2].at(podczeckShapes[2].size() - 1) << ','
                    << podczeckShapes[3].at(0) << ',' << podczeckShapes[3].at(podczeckShapes[3].size() * 0.25) << ',' << podczeckShapes[3].at(podczeckShapes[3].size() / 2) << ',' << podczeckShapes[3].at(podczeckShapes[3].size() * 0.75) << ',' << podczeckShapes[3].at(podczeckShapes[3].size() - 1) << ','
                    << podczeckShapes[4].at(0) << ',' << podczeckShapes[4].at(podczeckShapes[4].size() * 0.25) << ',' << podczeckShapes[4].at(podczeckShapes[4].size() / 2) << ',' << podczeckShapes[4].at(podczeckShapes[4].size() * 0.75) << ',' << podczeckShapes[4].at(podczeckShapes[4].size() - 1) << ',';
            }
            else {
                CSV << "0,0,"
                    << "0,0,0,0,0,"
                    << "0,0,0,0,0,"
                    << "0,0,0,0,0,"
                    << "0,0,0,0,0,"
                    << "0,0,0,0,0,"
                    << "0,0,0,0,0,"
                    << "0,0,0,0,0,";
            }
        }
    }
    catch (const std::exception& ex) {
        std::cout << "Exception: " << ex.what() << std::endl;
        exit(-1);
    }

    try {
        if (is_mean_curvature) {
            vcg::tri::UpdateFlags<TMesh>::FaceBorderFromFF(mesh);
            mesh.vert.EnableCurvature();
            mesh.vert.EnableQuality();
            if (vcg::tri::Clean<TMesh>::CountNonManifoldEdgeFF(mesh) > 0) {
                LOG << "[Wanring] Mesh has some not 2-manifold faces, Curvature computation requires manifoldness" << std::endl;
            }
            else {
                vcg::tri::UpdateCurvature<TMesh>::MeanAndGaussian(mesh);
                vcg::tri::UpdateQuality<TMesh>::VertexFromMeanCurvatureHG(mesh);
                vcg::Histogramf H;
                vcg::tri::Stat<TMesh>::ComputePerVertexQualityHistogram(mesh, H);
                H.SetRange(H.Percentile(0.1f), H.Percentile(0.9f), 100);
                std::stringstream _name;
                _name << filename << '_' << surface_name << run_timestamp << "_histogram.txt";
                H.FileWrite(_name.str());
                LOG << "-- Median curvature: " << H.Percentile(0.1f) << "," << H.Percentile(0.25f) << "," << H.Percentile(0.5f) << "," << H.Percentile(0.75f) << "," << H.Percentile(0.9f) << std::endl;
        else
            CSV << H.Percentile(0.1f) << "," << H.Percentile(0.25f) << "," << H.Percentile(0.5f) << "," << H.Percentile(0.75f) << "," << H.Percentile(0.9f) << ",";
            }
            mesh.vert.DisableCurvature();
            mesh.vert.DisableQuality();
        }
    }
    catch (const std::exception& ex) {
        std::cout << "Mean curvature exception: " << ex.what() << std::endl;
        exit(-1);
    }
    
    // Stage 2
    {
        // Report Volume and surface area
        int edgeNum = 0, edgeBorderNum = 0, edgeNonManifNum = 0;
        if (is_build_inverse) {
            vcg::tri::Clean<TMesh>::CountEdgeNum(inverse_mesh, edgeNum, edgeBorderNum, edgeNonManifNum);
            LOG << "[Pore connectivity]" << std::endl
                << "-- #Edge Border: " << edgeBorderNum << std::endl
                << "-- #Edge Non-manifold: " << edgeNonManifNum << std::endl
                << "-- Vertices: " << inverse_mesh.VN() << std::endl
                << "-- Faces: " << inverse_mesh.FN() << std::endl;
            bool watertight = (edgeBorderNum == 0) && (edgeNonManifNum == 0);
            bool pointcloud = (mesh.FN() == 0 && mesh.VN() != 0);
            if (!watertight || pointcloud) {
                LOG << "[Warning] Pore isn't conencted" << std::endl;
            }
        }

        vcg::tri::Clean<TMesh>::CountEdgeNum(mesh, edgeNum, edgeBorderNum, edgeNonManifNum);
        bool watertight = (edgeBorderNum == 0) && (edgeNonManifNum == 0);
        bool pointcloud = (mesh.FN() == 0 && mesh.VN() != 0);
        
        if (!pointcloud && watertight) {
            area2 = vcg::tri::Stat<TMesh>::ComputeMeshArea(mesh);
            volume2 = vcg::tri::Stat<TMesh>::ComputeMeshVolume(mesh);
            LOG << "[Scaffold properties]" << std::endl
                << "-- Volume: " << abs(volume2) << std::endl
                << "-- Surface Area: " << area2 << std::endl
                << "-- Porosity: " << 1 - abs(volume2 / volume1) << std::endl
                << "-- Surface Area ratio: " << area2 / area1 << std::endl;
            else
            CSV << abs(volume2) << ',' << area2 << ',' << 1 - abs(volume2 / volume1) << ',' << area2 / area1 << ','
                << mesh.VN() << ',' << mesh.FN() <<  std::endl;
        }
        else {
            LOG << "[Warning] The scaffold isn't a manifold. The grid_offset should be increased or enable --fix_self_intersect option" << std::endl;
            else
            CSV << ",,,," << mesh.VN() << ',' << mesh.FN() << std::endl;
        }

        if (!no_output) {
            LOG << "[Writing file] ";
            std::string filename2 = filename;
            filename.append("." + output_format);
            filename2.append("_inverse." + output_format);
            if (output_format == "ply") {
                vcg::tri::io::ExporterPLY<TMesh>::Save(mesh, filename.c_str(), false);
                if (is_build_inverse) vcg::tri::io::ExporterPLY<TMesh>::Save(inverse_mesh, filename2.c_str(), false);
            }
            else if (output_format == "obj") {
                vcg::tri::io::ExporterOBJ<TMesh>::Save(mesh, filename.c_str(), 0);
                if (is_build_inverse) vcg::tri::io::ExporterOBJ<TMesh>::Save(inverse_mesh, filename2.c_str(), 0);
            }
            else if (output_format == "off") {
                vcg::tri::io::ExporterOFF<TMesh>::Save(mesh, filename.c_str(), 0);
                if (is_build_inverse) vcg::tri::io::ExporterOFF<TMesh>::Save(inverse_mesh, filename2.c_str(), 0);
            }
            else if (output_format == "stl") {
                vcg::tri::io::ExporterSTL<TMesh>::Save(mesh, filename.c_str(), true, 0);
                if (is_build_inverse) vcg::tri::io::ExporterSTL<TMesh>::Save(inverse_mesh, filename2.c_str(), true, 0);
            }
            else if (output_format == "ctm") {
                vcg::tri::io::ExporterCTM<TMesh>::Save(mesh, filename.c_str(), 0, false, 1e-3);
                if (is_build_inverse) vcg::tri::io::ExporterCTM<TMesh>::Save(inverse_mesh, filename2.c_str(), 0, false, 1e-3);
            }
            LOG << "OK" << std::endl;
        }
        LOG << "[Finished]" << std::endl;
    }
    log_file.close();
    return 0;
}