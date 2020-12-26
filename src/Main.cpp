#include "Scaffolder_2.h"
#include "QuadricSimplification.h"

int main(int argc, char* argv[])
{
    // Define default parameters
    bool verbose = true;
    bool dirty = false;
    bool is_analysis_microstructure = false;
    bool is_export_microstructure = false;
    bool is_export_feret = false;
    bool is_build_inverse = false;
    bool is_fix_self_intersect = false;
    bool no_output = true;

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
    std::string filename = "";
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
            ("q,quiet", "Disable verbose output [default: false]")
            ("format", "Format of logging output [default: default]", cxxopts::value<std::string>(), "FORMAT (default, csv)")
            ("i,input", "Input file (STL)", cxxopts::value<std::string>(), "FILE")
            ("o,output", "Output filename with extension stl,ply,obj,off [default: out]", cxxopts::value<std::string>(), "FILENAME")
            ("output_inverse", "additional output inverse scaffold [default: false]")
            ("c,coff", "Angular frequency (pore size adjustment) default:PI", cxxopts::value<double>(), "DOUBLE")
            ("t,isolevel", "isolevel (porosity adjustment) [default: 0]", cxxopts::value<double>(), "DOUBLE")
            ("n,surface", "implicit surface: rectlinear, schwarzp, schwarzd, gyroid, double-p, double-d, double-gyroiod, lidinoid, schoen_iwp, neovius, bcc, tubular_g_ab, tubular_g_c [default: schwarzp]", cxxopts::value<std::string>(), "NAME")
            ("g,grid_size", "Grid size [default: 100]", cxxopts::value<size_t>(), "INT (0..60000)")
            ("params", "Combined parameters list", cxxopts::value<std::vector<std::string>>(), "surface[,coff,isolevel,grid_size,k_slice,k_polygon]")
            ("s,shell", "Outer thickness (layers) [default:0]", cxxopts::value<uint16_t>(), "INT (0..60000)")
            ("grid_offset", "[default:3]", cxxopts::value<uint16_t>(), "INT (0..60000)")
            ("m,microstructure", "Analysis microstructure with Slice contour technique ( [default: false]")
            ("export_microstructure", "Analysis microstructure and export the 2D contours (for debugging) [default: false]")
            ("k_slice", "K_slice: the number of slicing layers in each direction (used in microstructure analysis) [default: 100]", cxxopts::value<uint16_t>(), "INT (0..60000)")
            ("k_polygon", "K_polygon: the number of closest outer contour (used in microstructure analysis) [default: 4]", cxxopts::value<uint16_t>(), "INT (>0)")
            ("z,size_optimize", "Experimental Quadric simplification [default: 0]", cxxopts::value<double>(), "DOUBLE (0..1)")
            ("smooth_step", "Smooth with laplacian (default: 5)", cxxopts::value<uint16_t>(), "INT (0..60000)")
            ("dirty", "Disable autoclean [default false]")
            ("minimum_diameter", "used for removing small orphaned (between 0-1) [default: 0.25]", cxxopts::value<double>(), "DOUBLE (0..1)")
            ("fix_self_intersect", "Experimental fix self-intersect faces [default: false]");
        options.parse_positional({ "input", "output", "params" });
        bool isEmptyOption = (argc == 1);
        cxxopts::ParseResult result = options.parse(argc, argv);
        if (isEmptyOption || result.count("help")) {
            std::cout << options.help() << std::endl
                << "Example: " << std::endl
                << "  " << util::PathGetBasename(argv[0]) << " input.stl output.stl bcc,3.14159,0,100" << std::endl
                << "    " << "Generated BCC scaffold with w=3.14159 (PI), t=0, and grid size=100" << std::endl << std::endl;
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
            }
            if (output_format != "ply" && output_format != "obj" && output_format != "stl" && output_format != "off") {
                std::cout << "Invalid format: " << output_format << std::endl;
                return 1;
            }
            no_output = false;
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

        util::to_lower(surface_name);
        util::to_lower(output_format);

        if (surface_name == "rectlinear") {
            isolevel = 0;
        }
    }
    catch (const cxxopts::OptionException& ex) {
        std::cout << "Error parsing options: " << ex.what() << std::endl;
        return 1;
    }

    // Initialize output method via log file or standard output
    std::ofstream log_file;
    std::streambuf* log_buffer;
    if (verbose) {
        log_buffer = std::cout.rdbuf();
    }
    else {
        std::stringstream _name;
        _name << filename << '_' << surface_name << std::time(nullptr) << "." << ((log_format == SCAFFOLDER_FORMAT_CSV) ? "csv" : "txt");
        log_file.open(_name.str(), std::ofstream::out);
        log_buffer = log_file.rdbuf();
    }
    std::ostream log(log_buffer);
    std::ostream verbose_stream((verbose ? log_buffer : &util::null_buffer));

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
        << "-- Build Inverse: " << (is_build_inverse ? "True" : "False") << std::endl;
    else
    CSV << "Surface,coff,shell,thickness,grid_size,grid_offset,smooth_step,input_file,avg_min_feret,avg_max_feret,"
        << "min_min_feret,q1_min_feret,q2_min_feret,q3_min_feret,max_min_feret,"
        << "min_max_feret,q1_max_feret,q2_max_feret,q3_max_feret,max_max_feret,"
        << "min_square,q1_square,q2_square,q3_square,max_square,"
        << "min_circle,q1_circle,q2_circle,q3_circle,max_circle,"
        << "min_triangle,q1_triangle,q2_triangle,q3_triangle,max_triangle,"
        << "min_ellipse,q1_ellipse,q2_ellipse,q3_ellipse,max_ellipse,"
        << "min_elongation,q1_elongation,q2_elongation,q3_elongation,max_elongation,"
        << "volumn,surface_area,porosity,surface_area_ratio,vertices,faces" << std::endl
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
    {
        Eigen::MatrixXd V, Fxyz, IFxyz;
        Eigen::MatrixXi F;
        Eigen::RowVector3d V1min1;
        Eigen::RowVector3i grid_size;
        double grid_delta;
        try {
            // Read vertices and faces from STL 
            Eigen::MatrixXd V1;
            Eigen::MatrixXi F1;
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
            }

            {
                // Calculate the grid size parameters
                const Eigen::RowVector3d V1min = V1.colwise().minCoeff();
                const Eigen::RowVector3d V1max = V1.colwise().maxCoeff();
                const Eigen::RowVector3d L = V1max - V1min;
                const Eigen::RowVector3d delta = L / minimum_grid_size;
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
            
            Eigen::VectorXd W, D;
            {
                igl::FastWindingNumberBVH bvh;
                igl::fast_winding_number(V1, F1, 2, bvh);
                igl::fast_winding_number(bvh, 2, GV, W);
            }
            LOG << "OK" << std::endl;
            
            /*
            *  Shell generation
            *  Some application want to generate the porous scaffold with solid outer shell.
            *  this procedure utilize the winding number to thicken the outer shell of scaffold
            *  First, we find the border by search the neighbors around the inner points (MarkAndSweepNEighbor)
            *  If the some neighbors are outside, we mark this point as the border.
            *  Afterward, we use this border points to find the neighbors that are deeper inside and increase winding number
            *  of this neighbors by 1.
            *  The process is looped until shell level is matched by user-specific value
            *  Example:
            *  Initial 2D winding number grid:        After first loop:         After second loop:
            *       0 0 0 0 0 0 0 0 0 0             0 0 0 0 0 0 0 0 0 0         0 0 0 0 0 0 0 0 0 0
            *       0 1 1 1 1 1 0 1 0 0             0 3 3 3 3 3 0 3 0 0         0 3 3 3 3 3 0 3 0 0
            *       0 1 1 1 1 1 1 0 1 0             0 3 1 1 1 1 3 0 3 0         0 3 2 2 2 2 3 0 3 0
            *       0 0 1 1 1 1 0 0 1 0             0 0 3 1 1 3 0 0 3 0         0 0 3 2 2 3 0 0 3 0
            *       0 1 1 1 1 1 0 0 0 0             0 3 1 1 1 3 0 0 0 0         0 3 2 1 2 3 0 0 0 0
            *       0 0 1 1 1 0 0 0 0 0             0 0 3 1 3 0 0 0 0 0         0 0 3 2 3 0 0 0 0 0
            *       0 0 0 1 0 0 0 0 0 0             0 0 0 3 0 0 0 0 0 0         0 0 0 3 0 0 0 0 0 0
            */
            if (shell > 0) {
                LOG << "[Generate Shell] ";
                Queue_t queue;
                uint16_t shell_index = shell + 1;
                size_t grid_size_total = grid_size.prod();

                // find the border and set a new winding number value
                for (size_t index = 0; index < grid_size_total; index++) {
                    if (W(index) >= 0.8) {
                        index_type id;
                        indexToIJK(index, grid_size, id);
                        if (MarkAndSweepNeighbor(W, id, queue, grid_size, 0.5)) {
                            W(index) = shell_index;
                        }
                    }
                }
                shell_index--;
                // use the neighbors collected from previous step to find the deeper points and
                // set a new winding number value
                for (; shell_index > 1; shell_index--) {
                    Queue_t q(queue);
                    queue.clear();
                    for (Queue_t::iterator it = q.begin(); it != q.end(); ++it) {
                        // every neighbor
                        index_type id;
                        indexToIJK(it->first, grid_size, id);
                        if (MarkAndSweepNeighbor(W, id, queue, grid_size, shell_index, true)) {
                            W(it->first) = shell_index;
                        }
                    }
                }
                LOG << "OK" << std::endl;
            }
            
            LOG << "[Generating isosurface Fxyz] ";
            sol::state lua;
            Function* surface;
            sol::function f;
            // Initialize surface
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
                    lua.set("w", coff);
                    lua.set("t", isolevel);
                    lua.set_function("v", [](double x, double y, double z)
                    {
                            return 2*x;
                    });
                    lua_file();
                    f = lua["surface"];
                    surface = new LuaFunction(f);
                }
            }
            else {
               surface = isosurface(surface_name, isolevel);
            }
            Implicit_function fn(surface, coff);
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
            LOG << "OK" << std::endl;

            // Create scaffolder and inverse mesh by dual marching cube
            marching_cube(mesh, Fxyz, grid_size, V1min1, grid_delta, verbose, dirty);
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

            bool is_manifold = is_mesh_manifold(mesh);
            LOG << "-- is_manifold: " << is_manifold << std::endl;
            if (!is_manifold) {
                is_manifold = fix_non_manifold_vertices(mesh, minimum_diameter, 5, verbose_stream);
                is_manifold = is_manifold && fix_non_manifold_edges(mesh, minimum_diameter, 5, verbose_stream);
                LOG << "-- is_manifold: " << is_manifold << std::endl;
            }
            if (!is_manifold && !verbose) {
                LOG << "[Warning] Mesh is not two manifold" << std::endl;
            }
            if (log_format == SCAFFOLDER_FORMAT_DEFAULT) report_mesh(log, mesh);

            if (is_analysis_microstructure && is_manifold) {
                // Evaluating pore size by create 2D slice from 3D mesh
                std::string dir;
                std::stringstream filename;
                std::vector<double> minFeret;
                std::vector<double> maxFeret;
                std::vector<double> podczeckShapes[5];
                if (is_export_microstructure) {
                    // Initialize filename and directory for logging
                    // Get filename without extension from input_file.
                    // Ex., input.stl -> input
                    // then the folder `input_slice` will be created
                    dir = util::PathGetBasename(input_file) +"_slice";
                    util::make_dir(dir);
                }
                // init progress bar
                ProgressBar progress(grid_size.sum(), PROGRESS_BAR_COLUMN);
                
                // Start measure pore size by Slice contour technique
                vcg::tri::UpdateBounding<TMesh>::Box(mesh);
                vcg::Box3d bbox = mesh.bbox;
                vcg::Point3d dim = bbox.Dim();
                    
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
                    slice::Slice s = slice::incremental_slicing(mesh, k_slice, main_axis);
                    // update progress bar
                    progress += progress_major_step;
                    progress.display();
                    // Rodrigo's contour construct
                    slice::ContourSlice C = slice::contour_construct(s, main_axis);
                    // update progress bar
                    progress += progress_major_step;
                    progress.display();
                    // Measure pore size from C, and store the output in minFeret, maxFeret, and podczeckShapes
                    slice::measure_feret_and_shape(C, minFeret, maxFeret, podczeckShapes);
                    // update progress bar
                    progress += progress_major_step;
                    progress.display();
                    // If export flag was defined, convert the 2D contours into SVG
                    if (is_export_microstructure) {
                        filename.str(std::string());
                        unsigned int minor_step = C.size() / progress_major_step, step = minor_step;
                        // Loop foreach 2D contours
                        for (slice::ContourSlice::const_iterator cs = C.begin(); cs != C.end(); cs++) {
                            // Get array index from iterator
                            size_t index = (size_t)(cs - C.begin() + 1);
                            std::stringstream name(dir);
                            // Set the SVG filename to <axis>_<index>.svg
                            // Ex. x_0.svg, x_1.svg, ..., z_0.svg 
                            name << dir << "/" << axis << '_' << index << ".svg";
                            slice::write_svg(name.str(), *cs, dim[next_axis], dim[prev_axis], bbox.min[next_axis], bbox.min[prev_axis]);
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
                vcg::tri::io::ExporterSTL<TMesh>::Save(mesh, filename.c_str(), 0);
                if (is_build_inverse) vcg::tri::io::ExporterSTL<TMesh>::Save(inverse_mesh, filename2.c_str(), true, 0);
            }
            LOG << "OK" << std::endl;
        }
        LOG << "[Finished]" << std::endl;
        log_file.close();
    }
    return 0;
}