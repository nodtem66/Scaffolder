#include "PyScaffolder.hpp"

namespace py = pybind11;

PYBIND11_MODULE(PyScaffolder, m) {
    m.doc() = "PyScaffolder generate isosurface from implicit function";

    py::class_<PyScaffolder::PoreSize>(m, "PoreSize", py::dynamic_attr())
        .def(py::init<>())
        .def_readwrite("minFeret", &PyScaffolder::PoreSize::minFeret)
        .def_readwrite("maxFeret", &PyScaffolder::PoreSize::maxFeret);

    py::class_<PyScaffolder::MeshInfo>(m, "MeshInfo", py::dynamic_attr())
        .def(py::init<>())
        .def_readwrite("v", &PyScaffolder::MeshInfo::v)
        .def_readwrite("f", &PyScaffolder::MeshInfo::f)
        .def_readwrite("porosity", &PyScaffolder::MeshInfo::porosity)
        .def_readwrite("surface_area", &PyScaffolder::MeshInfo::surface_area)
        .def_readwrite("surface_area_ratio", &PyScaffolder::MeshInfo::surface_area_ratio);

    py::class_<PyScaffolder::Parameter>(m, "Parameter", py::dynamic_attr())
        .def(py::init<>())
        .def_readwrite("is_build_inverse", &PyScaffolder::Parameter::is_build_inverse)
        .def_readwrite("coff", &PyScaffolder::Parameter::coff)
        .def_readwrite("grid_offset", &PyScaffolder::Parameter::grid_offset)
        .def_readwrite("grid_size", &PyScaffolder::Parameter::grid_size)
        .def_readwrite("isolevel", &PyScaffolder::Parameter::isolevel)
        .def_readwrite("k_polygon", &PyScaffolder::Parameter::k_polygon)
        .def_readwrite("k_slice", &PyScaffolder::Parameter::k_slice)
        .def_readwrite("minimum_diameter", &PyScaffolder::Parameter::minimum_diameter)
        .def_readwrite("qsim_percent", &PyScaffolder::Parameter::qsim_percent)
        .def_readwrite("smooth_step", &PyScaffolder::Parameter::smooth_step)
        .def_readwrite("surface_name", &PyScaffolder::Parameter::surface_name);

    m.def("slice_test", &PyScaffolder::slice_test, "A function to slice input mesh into pore sizes",
        py::arg("vertices"),
        py::arg("faces"),
        py::arg("k_slice") = 100,
        py::arg("k_polygon") = 4,
        py::arg("direction") = 0,
        py::arg("callback") = py::none()
    );


    PyScaffolder::Parameter default_parameters;
    
    m.def("generate_mesh", &PyScaffolder::generate_mesh, "A function to generate isosurface from input mesh and parameters",
        py::arg("vertices"),
        py::arg("faces"),
        py::arg("params") = default_parameters,
        py::arg("callback") = py::none()
    );
}