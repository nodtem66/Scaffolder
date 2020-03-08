#include "Scaffolder_2.h"

int main(int argc, char* argv[])
{
    // Define default parameters
    bool verbose = true;
    bool dirty = false;
    bool is_analysis_microstructure = false;
    bool is_export_microstructure = false;
    bool is_export_feret = false;
    bool is_export_inverse = false;
    bool is_build_inverse = true;
    bool no_output = false;
    uint16_t grid_offset = 3;
    uint16_t shell = 0;
    uint16_t smooth_step = 5;
    uint8_t method = METHOD_SLICE_CONTOUR;
    double thickness = 0.0;
    size_t minimum_grid_size = 100;
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
            ("q,quiet", "Disable verbose output [default: false]")
            ("m,microstructure", "Analysis microstructure ( [default: false]")
            ("m1", "Export and analysis microstructure 1 (Image processing technique) [default: false]")
            ("m2", "Export and analysis microstructure 2 (Slice coutour technique) [default: false]")
            ("f,format", "Output format (OFF,PLY,STL,OBJ) [default: ply]", cxxopts::value<std::string>())
            ("i,input", "Input file (STL)", cxxopts::value<std::string>(), "FILE")
            ("o,output", "Output filename without extension [default: out]", cxxopts::value<std::string>(), "FILENAME")
            ("c,coff", "default:4*PI", cxxopts::value<double>(), "DOUBLE")
            ("s,shell", "[default:0]", cxxopts::value<uint16_t>(), "INT")
            ("n,surface", "rectlinear, schwarzp, schwarzd, gyroid, double-p, double-d, double-gyroiod, lidinoid, schoen_iwp, neovius, bcc, tubular_g_ab, tubular_g_c [default: schwarzp]", cxxopts::value<std::string>(), "NAME")
            ("t,thickness", "Thickness [default: 0]", cxxopts::value<double>(), "DOUBLE")
            ("g,grid_size", "Grid size [default: 100]", cxxopts::value<size_t>(), "INT")
            ("grid_offset", "[default:3]", cxxopts::value<uint16_t>(), "INT")
            ("smooth_step", "Smooth with laplacian (default: 5)", cxxopts::value<uint16_t>(), "INT")
            ("method", "Method of microstructure analysis: 0 (Image processing technique) or 1 (Slice contour technique) [default: 1]", cxxopts::value<uint8_t>(), "0,1")
            ("output_inverse", "additional output inverse scaffold")
            ("inverse", "Enable build inverse 3D scaffold (for pore connectivity analysis)")
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
        if (result.count("microstructure")) {
            is_analysis_microstructure = result["microstructure"].as<bool>();
            no_output = true;
        }
        if (result.count("method")) method = result["method"].as<uint8_t>();
        if (method < 0 && method > 1) {
            std::cout << "Invalid method option: " << method << ". The method option must be only 0 or 1, see --help";
            return 1;
        }
        if (result.count("m1")) {
            is_analysis_microstructure = result["m1"].as<bool>();
            is_export_microstructure = result["m1"].as<bool>();
            method = METHOD_IMAGE_PROCESS;
            no_output = true;
        }
        if (result.count("m2")) {
            is_analysis_microstructure = result["m2"].as<bool>();
            is_export_microstructure = result["m2"].as<bool>();
            method = METHOD_SLICE_CONTOUR;
            no_output = true;
        }
        if (result.count("output_inverse")) {
            is_export_inverse = result["output_inverse"].as<bool>();
        }
        if (result.count("output")) {
            filename = result["output"].as<std::string>();
            no_output = false;
        }
        if (result.count("format")) format = result["format"].as<std::string>();
        if (result.count("thickness")) thickness = result["thickness"].as<double>();
        if (result.count("grid_size")) minimum_grid_size = result["grid_size"].as<size_t>();
        if (result.count("grid_offset")) grid_offset = result["grid_offset"].as<uint16_t>();
        if (result.count("coff")) coff = result["coff"].as<double>();
        if (result.count("minimum_diameter")) minimum_diameter = result["minimum_diameter"].as<double>();
        if (result.count("surface")) surface = result["surface"].as<std::string>();
        if (result.count("shell")) shell = result["shell"].as<uint16_t>();
        if (result.count("smooth_step")) smooth_step = result["smooth_step"].as<uint16_t>();
        if (result.count("inverse")) is_build_inverse = result["inverse"].as<bool>();
        
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
    std::ofstream result;
    if (verbose) {
        std::cout << "[Scaffolder " << VERSION << "]" << std::endl
            << "-- Input file: " << input_file << std::endl
            << "-- Output file: " << filename << '.' << format << std::endl
            << "-- Surface: " << surface << std::endl
            << "-- Coff: " << coff << std::endl
            << "-- Thickness: " << thickness << std::endl
            << "-- Grid size: " << minimum_grid_size << std::endl
            << "--   Grid offset: " << grid_offset << std::endl
            << "--   Shell: " << shell << std::endl
            << "-- Autoclean: " << (dirty ? "False" : "True") << std::endl
            << "--   Minimum diameter: " << 100 * minimum_diameter << "%" << std::endl
            << "--   Smooth step: " << smooth_step << std::endl
            << "-- Analysis microstructure: " << (is_analysis_microstructure ? "True" : "False") << std::endl
            << "--   Method :" << (method == METHOD_IMAGE_PROCESS ? "Image Processing" : "Contour Slicing") << std::endl
            << "--   Export microstructure: " << (is_analysis_microstructure ? "True" : "False") << std::endl
            << "-- Build Inverse: " << (is_build_inverse ? "True" : "False") << std::endl
            << "--   Export Inverse: " << (is_export_inverse ? "True" : "False") << std::endl;
    }
    else {
        std::stringstream _name;
        _name << filename << '_' << surface << std::time(nullptr) << ".txt";
        result.open(_name.str(), std::ofstream::out);
        // Print header
        result << "Surface,coff,shell,thickness,grid_size,grid_offset,smooth_step,input_file,avg_min_feret,avg_max_feret,"
            << "min_min_feret,q1_min_feret,q2_min_feret,q3_min_feret,max_min_feret,"
            << "min_max_feret,q1_max_feret,q2_max_feret,q3_max_feret,max_max_feret,"
            << "vol,surface_area,porosity,surface_area_ratio,#vertices,#faces,"
            << "min_square,q1_square,q2_square,q3_square,max_square,"
            << "min_circle,q1_circle,q2_circle,q3_circle,max_circle,"
            << "min_triangle,q1_triangle,q2_triangle,q3_triangle,max_triangle,"
            << "min_ellipse,q1_ellipse,q2_ellipse,q3_ellipse,max_ellipse,"
            << "min_elongation,q1_elongation,q2_elongation,q3_elongation,max_elongation,"
            << "volumn,surface_area,porosity,surface_area_ratio,vertices,faces" << std::endl;
        result << surface << ',' << coff << ',' << shell << ',' << thickness << ',' << minimum_grid_size << ',' << grid_offset << ',' << smooth_step << ','
            << input_file << ',';
    }

    // Stage 1:
    TMesh mesh, inverse_mesh;
    uint64_t starttime, endtime;
    double volume1 = eps, volume2 = eps;
    double area1 = eps, area2 = eps;
    {
        Eigen::MatrixXd V, Fxyz, IFxyz;
        Eigen::MatrixXi F;
        Eigen::RowVector3d V1min1;
        Eigen::RowVector3i grid_size;
        double grid_delta;
        try {
            // Read FROM STL
            Eigen::MatrixXd V1;
            Eigen::MatrixXi F1;
            {
                TMesh stl;
                int loadmark = 0;
                vcg::tri::io::ImporterSTL<TMesh>::Open(stl, input_file.c_str(), loadmark);
                vcg::tri::Clean<TMesh>::RemoveDuplicateFace(stl);
                vcg::tri::Clean<TMesh>::RemoveDuplicateVertex(stl);
                vcg::tri::Clean<TMesh>::RemoveUnreferencedVertex(stl);
                vcg::tri::UpdateTopology<TMesh>::FaceFace(stl);
                vcg::tri::Clean<TMesh>::RemoveZeroAreaFace(stl);
                vcg::tri::Clean<TMesh>::RemoveNonManifoldFace(stl);
                vcg::tri::Clean<TMesh>::RemoveUnreferencedVertex(stl);
                vcg::tri::UpdateTopology<TMesh>::FaceFace(stl);
                // Report Volume and surface area
                int edgeNum = 0, edgeBorderNum = 0, edgeNonManifNum = 0;
                vcg::tri::Clean<TMesh>::CountEdgeNum(stl, edgeNum, edgeBorderNum, edgeNonManifNum);
                bool watertight = (edgeBorderNum == 0) && (edgeNonManifNum == 0);
                bool pointcloud = (stl.fn == 0 && stl.vn != 0);
                if (!pointcloud && watertight) {
                    area1 = vcg::tri::Stat<TMesh>::ComputeMeshArea(stl);
                    volume1 = vcg::tri::Stat<TMesh>::ComputeMeshVolume(stl);
                }
                else {
                    std::cout << "Error: Input file is not valid" << std::endl
                        << "-- Watertight: " << watertight << (watertight ? "[Valid]" : "[Invalid]") << std::endl
                        << "-- Point cloud: " << pointcloud << (!pointcloud ? "[Valid]" : "[Invalid]") << std::endl;
                    return 1;
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
                if (verbose)
                    std::cout << "-- Bounding Box: "
                    << V1min.format(CleanFmt) << ' ' << V1max.format(CleanFmt) << std::endl
                    << "-- Length: " << L.format(CleanFmt) << std::endl
                    << "-- Grid delta: " << grid_delta << std::endl;
            }
            Implicit_function fn(isosurface(surface, thickness), coff);
            if (verbose) std::cout << "[Generating grid] ";
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
            if (verbose) 
                std::cout
                    << "OK" << std::endl
                    << "-- Grid size: " << grid_size.prod() << " " << grid_size.format(CleanFmt) << std::endl;

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
                size_t grid_size_total = grid_size.prod();
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
                if (verbose) std::cout << "OK" << std::endl;
            }
            
            if (verbose) std::cout << "[Generating isosurface Fxyz] ";
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
                        // Winding field
                        if (w < 0.8) {// Outside
                            Fxyz(index) = 1.0;
                            IFxyz(index) = 1.0;
                        }
                        else if (w >= 1.1) {// Shell
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
            if (verbose) std::cout << "OK" << std::endl;

            // Create scaffolder and inverse mesh by dual marching cube
            marching_cube(mesh, Fxyz, grid_size, V1min1, grid_delta, verbose, dirty);
            if (is_build_inverse) marching_cube(inverse_mesh, IFxyz, grid_size, V1min1, grid_delta, false, false);

            if (!dirty) {
                clean_mesh(mesh, minimum_diameter, smooth_step, verbose);
                if (is_build_inverse) clean_mesh(inverse_mesh, minimum_diameter, 0, false);
            }
            bool is_manifold = is_mesh_manifold(mesh);
            if (!is_manifold && !verbose) {
                std::cout << "[Warning] Mesh is not two manifold" << std::endl;
            }
            if (verbose) report_mesh(mesh);

            if (is_analysis_microstructure && is_manifold) {
                // Evaluating pore size by create 2D slice 8-bit image (0-blacks're pores and 255-whites're grains)
                // Then an 8-bit image become a binary image by image thresholding (value of 150)
                // The binary imaege'll be labeled and finally evaluated the feret diameter by chain coding
            
                // init filename
                std::string dir;
                std::stringstream filename;
                std::vector<dip::dfloat> minFeret;
                std::vector<dip::dfloat> maxFeret;
                std::vector<dip::dfloat> podczeckShapes[5];
                if (is_export_microstructure) {
                    size_t firstindex = input_file.find_last_of("/\\");
                    firstindex = firstindex == string::npos ? 0 : firstindex + 1;
                    size_t lastindex = input_file.find_last_of(".");
                    dir = input_file.substr(firstindex, lastindex - firstindex) + "_slice";
                    make_dir(dir);
                }
                // init progress bar
                ProgressBar progress(grid_size.sum(), PROGRESS_BAR_COLUMN);
                
                if (method == METHOD_IMAGE_PROCESS) {
                    // init DIP MeasurementTool and DIP image
                    dip::MeasurementTool tool;
                    dip::Image img2d({ (dip::uint64)grid_size.maxCoeff(), (dip::uint64)grid_size.maxCoeff() }, 1, dip::DT_UINT8);
                    dip::uint8* data = static_cast<dip::uint8*>(img2d.Origin());
                    char axis = 'x';
                    for (uint8_t main_axis = 0; main_axis < 3; main_axis++) {
                        // main_axis: 0 = x, 1 = y, 2 = z
                        size_t index;
                        uint8_t axis2 = (main_axis + 1) % 3;
                        uint8_t axis3 = (main_axis + 2) % 3;
                        for (size_t k = grid_offset; k < grid_size(main_axis) - (size_t)(grid_offset); k++) {
                            img2d.Fill(200);
                            for (size_t j = grid_offset; j < grid_size(axis2) - (size_t)(grid_offset); j++) {
                                for (size_t i = grid_offset; i < grid_size(axis3) - (size_t)(grid_offset); i++) {
                                    if (main_axis == 0)
                                        index = k + grid_size(0) * (j + grid_size(1) * i);
                                    else if (main_axis == 1)
                                        index = i + grid_size(0) * (k + grid_size(1) * j);
                                    else
                                        index = j + grid_size(0) * (i + grid_size(1) * k);
                                    if (W(index) >= 0.8) {// inside STL mesh
                                        if (Fxyz(index) > eps2)
                                            data[i + grid_size(axis3) * j] = 0;
                                    }
                                }
                            }
                            // Measurement Feret diameter
                            dip::Image label = dip::Label(img2d < 50, 2);
                            dip::Measurement msr = tool.Measure(label, img2d, { "Feret", "PodczeckShapes" }, {}, 2);
                            dip::Measurement::IteratorFeature it = msr["Feret"];
                            dip::Measurement::IteratorFeature::Iterator feret = it.FirstObject();
                            while (feret) {
                                // From ref: https://diplib.github.io/diplib-docs/features.html#size_features_Feret
                                maxFeret.push_back(feret[2] * grid_delta);
                                minFeret.push_back(feret[1] * grid_delta);
                                ++feret;
                            }
                            it = msr["PodczeckShapes"];
                            dip::Measurement::IteratorFeature::Iterator shape = it.FirstObject();
                            while (shape) {
                                podczeckShapes[0].push_back(shape[0]);
                                podczeckShapes[1].push_back(shape[1]);
                                podczeckShapes[2].push_back(shape[2]);
                                podczeckShapes[3].push_back(shape[3]);
                                podczeckShapes[4].push_back(shape[4]);
                                ++shape;
                            }
                            if (is_export_microstructure) {
                                filename.str(std::string());
                                filename << dir << '/' << axis << '_' << k << ".jpg";
                                dip::ImageWriteJPEG(img2d, filename.str());
                            }
                            ++progress;
                            progress.display();
                        }
                        axis++;
                    }
                } // End Image processing
                else if (method == METHOD_SLICE_CONTOUR) {
                    vcg::tri::UpdateBounding<TMesh>::Box(mesh);
                    vcg::Box3d bbox = mesh.bbox;
                    vcg::Point3d dim = bbox.Dim();
                    
                    char axis = 'x';
                    progress.display();
                    for (uint8_t main_axis = 0; main_axis < 3; main_axis++) {
                        // main_axis: 0 = x, 1 = y, 2 = z
                        size_t index;
                        int next_axis = ((main_axis + 1) % 3);
                        int prev_axis = ((main_axis + 2) % 3);
                        unsigned int progress_major_step = grid_size[main_axis] / 4;
                        slice::Slice s = slice::incremental_slicing(mesh, grid_size(main_axis), main_axis);
                        progress += progress_major_step;
                        progress.display();
                        slice::ContourSlice C = slice::contour_construct(s, main_axis);
                        progress += progress_major_step;
                        progress.display();
                        slice::measure_feret_and_shape(C, minFeret, maxFeret, podczeckShapes);
                        progress += progress_major_step;
                        progress.display();
                        if (is_export_microstructure) {
                            filename.str(std::string());
                            unsigned int minor_step = C.size() / progress_major_step, step = minor_step;
                            for (slice::ContourSlice::const_iterator cs = C.begin(); cs != C.end(); cs++) {
                                size_t index = (size_t)(cs - C.begin() + 1);
                                std::stringstream name(dir);
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
                } // End Slice contour
                progress.done();
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
                    if (verbose) {
                        std::cout << "[Microstructure] " << std::endl
                            << "-- Avg Min Feret: " << std::accumulate(minFeret.begin(), minFeret.end(), 0.0) / minFeret.size() << std::endl
                            << "-- Avg Max Feret: " << std::accumulate(maxFeret.begin(), maxFeret.end(), 0.0) / maxFeret.size() << std::endl
                            << "-- Min Feret: [" << minFeret.at(0) << " " << minFeret.at(minFeret.size() * 0.25) << " " << minFeret.at(minFeret.size() / 2) << " " << minFeret.at(minFeret.size() * 0.75) << " " << minFeret.at(minFeret.size() - 1) << "]" << std::endl
                            << "-- Max Feret: [" << maxFeret.at(0) << " " << maxFeret.at(maxFeret.size() * 0.25) << " " << maxFeret.at(maxFeret.size() / 2) << " " << maxFeret.at(maxFeret.size() * 0.75) << " " << maxFeret.at(maxFeret.size() - 1) << "]" << std::endl
                            << "-- Square Similarity: [" << podczeckShapes[0].at(0) << " " << podczeckShapes[0].at(podczeckShapes[0].size() * 0.25) << " " << podczeckShapes[0].at(podczeckShapes[0].size() / 2) << " " << podczeckShapes[0].at(podczeckShapes[0].size() * 0.75) << " " << podczeckShapes[0].at(podczeckShapes[0].size() - 1) << "]" << std::endl
                            << "-- Circle Similarity: [" << podczeckShapes[1].at(0) << " " << podczeckShapes[1].at(podczeckShapes[1].size() * 0.25) << " " << podczeckShapes[1].at(podczeckShapes[1].size() / 2) << " " << podczeckShapes[1].at(podczeckShapes[1].size() * 0.75) << " " << podczeckShapes[1].at(podczeckShapes[1].size() - 1) << "]" << std::endl
                            << "-- Triangle Similarity: [" << podczeckShapes[2].at(0) << " " << podczeckShapes[2].at(podczeckShapes[2].size() * 0.25) << " " << podczeckShapes[2].at(podczeckShapes[2].size() / 2) << " " << podczeckShapes[2].at(podczeckShapes[2].size() * 0.75) << " " << podczeckShapes[2].at(podczeckShapes[2].size() - 1) << "]" << std::endl
                            << "-- Ellipse Similarity: [" << podczeckShapes[3].at(0) << " " << podczeckShapes[3].at(podczeckShapes[3].size() * 0.25) << " " << podczeckShapes[3].at(podczeckShapes[3].size() / 2) << " " << podczeckShapes[3].at(podczeckShapes[3].size() * 0.75) << " " << podczeckShapes[3].at(podczeckShapes[3].size() - 1) << "]" << std::endl
                            << "-- Elongation Similarity: [" << podczeckShapes[4].at(0) << " " << podczeckShapes[4].at(podczeckShapes[4].size() * 0.25) << " " << podczeckShapes[4].at(podczeckShapes[4].size() / 2) << " " << podczeckShapes[4].at(podczeckShapes[4].size() * 0.75) << " " << podczeckShapes[4].at(podczeckShapes[4].size() - 1) << "]" << std::endl;
                    }
                    else {
                        result << std::accumulate(minFeret.begin(), minFeret.end(), 0.0) / minFeret.size() << ','
                            << std::accumulate(maxFeret.begin(), maxFeret.end(), 0.0) / maxFeret.size() << ','
                            << minFeret.at(0) << ',' << minFeret.at(minFeret.size() * 0.25) << ',' << minFeret.at(minFeret.size() / 2) << ',' << minFeret.at(minFeret.size() * 0.75) << ',' << minFeret.at(minFeret.size() - 1) << ','
                            << maxFeret.at(0) << ',' << maxFeret.at(maxFeret.size() * 0.25) << ',' << maxFeret.at(maxFeret.size() / 2) << ',' << maxFeret.at(maxFeret.size() * 0.75) << ',' << maxFeret.at(maxFeret.size() - 1) << ','
                            << podczeckShapes[0].at(0) << ',' << podczeckShapes[0].at(podczeckShapes[0].size() * 0.25) << ',' << podczeckShapes[0].at(podczeckShapes[0].size() / 2) << ',' << podczeckShapes[0].at(podczeckShapes[0].size() * 0.75) << ',' << podczeckShapes[0].at(podczeckShapes[0].size() - 1) << ','
                            << podczeckShapes[1].at(0) << ',' << podczeckShapes[1].at(podczeckShapes[1].size() * 0.25) << ',' << podczeckShapes[1].at(podczeckShapes[1].size() / 2) << ',' << podczeckShapes[1].at(podczeckShapes[1].size() * 0.75) << ',' << podczeckShapes[1].at(podczeckShapes[1].size() - 1) << ','
                            << podczeckShapes[2].at(0) << ',' << podczeckShapes[2].at(podczeckShapes[2].size() * 0.25) << ',' << podczeckShapes[2].at(podczeckShapes[2].size() / 2) << ',' << podczeckShapes[2].at(podczeckShapes[2].size() * 0.75) << ',' << podczeckShapes[2].at(podczeckShapes[2].size() - 1) << ','
                            << podczeckShapes[3].at(0) << ',' << podczeckShapes[3].at(podczeckShapes[3].size() * 0.25) << ',' << podczeckShapes[3].at(podczeckShapes[3].size() / 2) << ',' << podczeckShapes[3].at(podczeckShapes[3].size() * 0.75) << ',' << podczeckShapes[3].at(podczeckShapes[3].size() - 1) << ','
                            << podczeckShapes[4].at(0) << ',' << podczeckShapes[4].at(podczeckShapes[4].size() * 0.25) << ',' << podczeckShapes[4].at(podczeckShapes[4].size() / 2) << ',' << podczeckShapes[4].at(podczeckShapes[4].size() * 0.75) << ',' << podczeckShapes[4].at(podczeckShapes[4].size() - 1) << ',';
                    }
                }
                else {
                    result << "0,0,"
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
            std::cout << strerror(errno) << endl;
            return errno;
        }
    }

    // Stage 2: Result
    {
        // Report Volume and surface area
        int edgeNum = 0, edgeBorderNum = 0, edgeNonManifNum = 0;
        if (is_build_inverse) {
            vcg::tri::Clean<TMesh>::CountEdgeNum(inverse_mesh, edgeNum, edgeBorderNum, edgeNonManifNum);
            if (verbose)
                std::cout << "[Pore connectivity]" << std::endl
                << "-- #Edge Border: " << edgeBorderNum << std::endl
                << "-- #Edge Non-manifold: " << edgeNonManifNum << std::endl
                << "-- Vertices: " << inverse_mesh.VN() << std::endl
                << "-- Faces: " << inverse_mesh.FN() << std::endl;
            bool watertight = (edgeBorderNum == 0) && (edgeNonManifNum == 0);
            bool pointcloud = (mesh.fn == 0 && mesh.vn != 0);
            if (!watertight || pointcloud) std::cout << "[Warning] Pore isn't conencted" << std::endl;
        }

        vcg::tri::Clean<TMesh>::CountEdgeNum(mesh, edgeNum, edgeBorderNum, edgeNonManifNum);
        bool watertight = (edgeBorderNum == 0) && (edgeNonManifNum == 0);
        bool pointcloud = (mesh.fn == 0 && mesh.vn != 0);
        
        if (!pointcloud && watertight) {
            area2 = vcg::tri::Stat<TMesh>::ComputeMeshArea(mesh);
            volume2 = vcg::tri::Stat<TMesh>::ComputeMeshVolume(mesh);
            if (verbose)
                std::cout << "[Scaffold properties]" << std::endl
                << "-- Volume: " << abs(volume2) << std::endl
                << "-- Surface Area: " << area2 << std::endl
                << "-- Porosity: " << abs(volume2 / volume1) << std::endl
                << "-- Surface Area ratio: " << area2 / area1 << std::endl;
            else
                result
                << abs(volume2) << ',' << area2 << ',' << abs(volume2 / volume1) << ',' << area2 / area1 << ','
                << mesh.VN() << ',' << mesh.FN() <<  std::endl;
        }
        else {
            if (!verbose) {
                result << ",,,," << mesh.VN() << ',' << mesh.FN() << std::endl;
            }
            std::cout << "[Warning] The scaffolder isn't a manifold. The grid_offset should have been increased" << std::endl;
        }
        
        if (!no_output) {
            if (verbose) std::cout << "[Writing file] ";
            std::string filename2 = filename;
            filename.append("." + format);
            filename2.append("_inverse." + format);
            if (format == "ply") {
                vcg::tri::io::ExporterPLY<TMesh>::Save(mesh, filename.c_str(), false);
                if (is_build_inverse && is_export_inverse) vcg::tri::io::ExporterPLY<TMesh>::Save(inverse_mesh, filename2.c_str(), false);
            }
            else if (format == "obj") {
                vcg::tri::io::ExporterOBJ<TMesh>::Save(mesh, filename.c_str(), 0);
                if (is_build_inverse && is_export_inverse) vcg::tri::io::ExporterOBJ<TMesh>::Save(inverse_mesh, filename2.c_str(), 0);
            }
            else if (format == "off") {
                vcg::tri::io::ExporterOFF<TMesh>::Save(mesh, filename.c_str(), 0);
                if (is_build_inverse && is_export_inverse) vcg::tri::io::ExporterOFF<TMesh>::Save(inverse_mesh, filename2.c_str(), 0);
            }
            else if (format == "stl") {
                vcg::tri::io::ExporterSTL<TMesh>::Save(mesh, filename.c_str(), 0);
                if (is_build_inverse && is_export_inverse) vcg::tri::io::ExporterSTL<TMesh>::Save(inverse_mesh, filename2.c_str(), 0);
            }
            if (verbose) std::cout << "OK" << std::endl;
        }
        if (verbose) std::cout << "[Finished]" << std::endl;
        result.close();
    }
    return 0;
}