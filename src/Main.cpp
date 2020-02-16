#include "Scaffolder_2.h"

int main(int argc, char* argv[])
{
    // Define default parameters
    bool verbose = true;
    bool dirty = false;
    bool is_analysis_microstructure = false;
    bool is_export_microstructure = false;
    bool is_export_feret = false;
    bool is_export_hdf5 = false;
    bool no_output = false;
    uint16_t grid_offset = 2;
    uint16_t shell = 0;
    uint16_t smooth_step = 5;
    double thickness = 0.5;
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
            ("m,microstructure", "Analysis microstructure [default: false]")
            ("m2", "Export and analysis microstructure [default: false]")
            ("m3", "Export feret data")
            ("f,format", "Output format (OFF,PLY,STL,OBJ) [default: ply]", cxxopts::value<std::string>())
            ("i,input", "Input file (STL)", cxxopts::value<std::string>(), "FILE")
            ("o,output", "Output filename without extension [default: out]", cxxopts::value<std::string>(), "FILENAME")
            ("c,coff", "default:4*PI", cxxopts::value<double>(), "DOUBLE")
            ("s,shell", "[default:0]", cxxopts::value<uint16_t>(), "INT")
            ("n,surface", "rectlinear, schwarzp, schwarzd, gyroid, double-p, double-d, double-gyroiod, lidinoid, schoen_iwp, neovius, bcc, tubular_g_ab, tubular_g_c [default: schwarzp]", cxxopts::value<std::string>(), "NAME")
            ("t,thickness", "Thickness [default: 1.0]", cxxopts::value<double>(), "DOUBLE")
            ("g,grid_size", "Grid size [default: 100]", cxxopts::value<size_t>(), "INT")
            ("grid_offset", "[default:2]", cxxopts::value<uint16_t>(), "INT")
            ("smooth_step", "Smooth with laplacian (default: 5)", cxxopts::value<uint16_t>(), "INT")
            ("hdf5", "export the result in hdf5 format [default: false]")
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
        if (result.count("m2")) {
            is_analysis_microstructure = result["m2"].as<bool>();
            is_export_microstructure = result["m2"].as<bool>();
            no_output = true;
        }
        if (result.count("m3")) {
            is_export_feret = result["m3"].as<bool>();
            is_analysis_microstructure = is_export_feret;
        }
        if (result.count("hdf5")) {
            is_export_hdf5 = result["hdf5"].as<bool>();
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
            << "-- Grid size: " << minimum_grid_size << std::endl
            << "-- Grid offset: " << grid_offset << std::endl
            << "-- Smooth step: " << smooth_step << std::endl
            << "-- Autoclean: " << (dirty ? "False" : "True") << std::endl
            << "-- Minimum diameter: " << 100 * minimum_diameter << "%" << std::endl
            << "-- Analysis microstructure: " << (is_analysis_microstructure || is_export_microstructure ? "True" : "False") << std::endl
            << "-- Export HDF5: " << (is_export_hdf5 ? "True" : "False") << std::endl;
    }

    // Stage 1:
    TMesh mesh, inverse_mesh;
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
                if (!dirty) {
                    vcg::tri::Clean<TMesh>::RemoveDuplicateFace(stl);
                    vcg::tri::Clean<TMesh>::RemoveDuplicateVertex(stl);
                    vcg::tri::Clean<TMesh>::RemoveUnreferencedVertex(stl);
                    vcg::tri::UpdateTopology<TMesh>::FaceFace(stl);
                    vcg::tri::Clean<TMesh>::RemoveZeroAreaFace(stl);
                    vcg::tri::Clean<TMesh>::RemoveNonManifoldFace(stl);
                    vcg::tri::Clean<TMesh>::RemoveUnreferencedVertex(stl);
                    vcg::tri::UpdateTopology<TMesh>::FaceFace(stl);
                }
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
                // Convert STL VCG::TMesh
                // Vertices to Eigen matrixXd V1
                V1.resize(stl.VN(), 3);
                size_t i = 0;
                std::vector<size_t> vertexId(stl.vert.size());
                for (TMesh::VertexIterator it = stl.vert.begin(); it != stl.vert.end(); ++it) if (!it->IsD()) {
                    vertexId[it - stl.vert.begin()] = i;
                    vcg::Point3d point = it->P();
                    V1(i, 0) = point[0];
                    V1(i, 1) = point[1];
                    V1(i, 2) = point[2];
                    i++;
                }
                // Faces to Eigen matrixXi F1
                i = 0;
                F1.resize(stl.FN(), stl.face.begin()->VN());
                for (TMesh::FaceIterator it = stl.face.begin(); it != stl.face.end(); ++it) if (!it->IsD()) {
                    for (int k = 0; k < it->VN(); k++) {
                        F1(i, k) = vertexId[vcg::tri::Index(stl, it->V(k))];
                    }
                    i++;
                }
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
                if (verbose) std::cout << "-- Bounding Box: " << V1min.format(CleanFmt) << ' ' << V1max.format(CleanFmt) << std::endl;
                if (verbose) std::cout << "-- Length: " << L << std::endl;
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
                    << "-- Grid size: " << grid_size.prod() << grid_size.format(CleanFmt) << std::endl;

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

            // Evaluating pore size by create 2D slice 8-bit image (0-blacks're pores and 255-whites're grains)
            // Then an 8-bit image become a binary image by image thresholding (value of 150)
            // The binary imaege'll be labeled and finally evaluated the feret diameter by chain coding
            if (is_analysis_microstructure) {
                // init filename
                std::string dir;
                std::stringstream filename;
                std::vector<dip::dfloat> minFeret;
                std::vector<dip::dfloat> maxFeret;
                if (is_export_microstructure) {
                    size_t firstindex = input_file.find_last_of("/\\");
                    firstindex = firstindex == string::npos ? 0 : firstindex + 1;
                    size_t lastindex = input_file.find_last_of(".");
                    dir = input_file.substr(firstindex, lastindex - firstindex) + "_slice";
                    make_dir(dir);
                }
                // init progress bar
                ProgressBar progress(grid_size.sum(), PROGRESS_BAR_COLUMN);
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
                    for (size_t k = 0; k < grid_size(main_axis); k++) {
                        img2d.Fill(255);
                        for (size_t j = 0; j < grid_size(axis2); j++) {
                            for (size_t i = 0; i < grid_size(axis3); i++) {
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
                        dip::Image label = dip::Label(img2d < 150, 2);
                        dip::Measurement msr = tool.Measure(label, img2d, { "Feret" }, {}, 2);
                        dip::Measurement::IteratorFeature it = msr["Feret"];
                        dip::Measurement::IteratorFeature::Iterator feret = it.FirstObject();
                        while (feret) {
                            maxFeret.push_back(feret[0]*grid_delta);
                            minFeret.push_back(feret[1]*grid_delta);
                            ++feret;
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
                if (is_export_feret) {
                    std::string f = filename.str() + ".csv";
                    std::ofstream fs(f.c_str(), std::ofstream::out);
                    fs << "MIN";
                    for (std::vector<dip::dfloat>::iterator it = minFeret.begin(); it != minFeret.end(); it++) {
                        fs << "," << *it;
                    }
                    fs << std::endl;
                    fs << "MAX";
                    for (std::vector<dip::dfloat>::iterator it = maxFeret.begin(); it != maxFeret.end(); it++) {
                        fs << "," << *it;
                    }
                    fs << std::endl;
                    fs.close();
                }
                progress.done();
                if (verbose) {
                    std::cout << "[Microstructure] " << std::endl
                        << "-- Avg Min Feret: " << std::accumulate(minFeret.begin(), minFeret.end(), 0.0) / minFeret.size() << std::endl
                        << "-- Avg Max Feret: " << std::accumulate(maxFeret.begin(), maxFeret.end(), 0.0) / maxFeret.size() << std::endl;
                }
            }
        }
        catch (const std::exception& ex) {
            std::cout << "Exception: " << ex.what() << std::endl;
            std::cout << strerror(errno) << endl;
            return errno;
        }

        marching_cube(mesh, Fxyz, grid_size, V1min1, grid_delta, verbose, dirty);
        marching_cube(inverse_mesh, IFxyz, grid_size, V1min1, grid_delta, false, false);
    }

    // Stage 2: Cleaning
    {
        if (!dirty) {
            clean_mesh(mesh, minimum_diameter, smooth_step, verbose);
            clean_mesh(inverse_mesh, minimum_diameter, smooth_step, false);
        }

        // Report Volume and surface area
        int edgeNum = 0, edgeBorderNum = 0, edgeNonManifNum = 0;
        vcg::tri::Clean<TMesh>::CountEdgeNum(inverse_mesh, edgeNum, edgeBorderNum, edgeNonManifNum);

        if (verbose)
            std::cout << "[Pore connectivity]" << std::endl
            << "-- #Edge Border: " << edgeBorderNum << std::endl
            << "-- #Edge Non-manifold: " << edgeNonManifNum << std::endl;

        vcg::tri::Clean<TMesh>::CountEdgeNum(mesh, edgeNum, edgeBorderNum, edgeNonManifNum);
        bool watertight = (edgeBorderNum == 0) && (edgeNonManifNum == 0);
        bool pointcloud = (mesh.fn == 0 && mesh.vn != 0);
        
        if (!pointcloud && watertight) {
            area2 = vcg::tri::Stat<TMesh>::ComputeMeshArea(mesh);
            volume2 = vcg::tri::Stat<TMesh>::ComputeMeshVolume(mesh);
            if (verbose)
                std::cout 
                << "-- Volume: " << volume2 << std::endl
                << "-- Surface Area: " << area2 << std::endl
                << "-- Porosity: " << volume2 / volume1 << std::endl
                << "-- Surface Area ratio: " << area2 / area1 << std::endl;
        }
        
        if (!no_output) {
            if (verbose) std::cout << "[Writing file] ";
            std::string filename2 = filename;
            filename.append("." + format);
            filename2.append("_inverse." + format);
            if (format == "ply") {
                vcg::tri::io::ExporterPLY<TMesh>::Save(mesh, filename.c_str(), false);
                vcg::tri::io::ExporterPLY<TMesh>::Save(inverse_mesh, filename2.c_str(), false);
            }
            else if (format == "obj") {
                vcg::tri::io::ExporterOBJ<TMesh>::Save(mesh, filename.c_str(), 0);
                vcg::tri::io::ExporterOBJ<TMesh>::Save(inverse_mesh, filename2.c_str(), 0);
            }
            else if (format == "off") {
                vcg::tri::io::ExporterOFF<TMesh>::Save(mesh, filename.c_str(), 0);
                vcg::tri::io::ExporterOFF<TMesh>::Save(inverse_mesh, filename2.c_str(), 0);
            }
            else if (format == "stl") {
                vcg::tri::io::ExporterSTL<TMesh>::Save(mesh, filename.c_str(), 0);
                vcg::tri::io::ExporterSTL<TMesh>::Save(inverse_mesh, filename2.c_str(), 0);
            }
        }
        if (verbose) std::cout << "OK" << std::endl << "[Finished]" << std::endl;
    }
    return 0;
}