// Copyright 1990-2026 by Thomas M. Breuel
// Licensed under the Apache License, Version 2.0 (see LICENSE)

#include <pybind11/pybind11.h>

#include "cedges.h"
#include "rast.h"

namespace py = pybind11;

PYBIND11_MODULE(rast, m) {
  m.doc() = "RAST family of branch-and-bound feature matching algorithms";

  py::class_<LinesP2D>(m, "LinesP2D")
      .def("set_maxresults", &LinesP2D::set_maxresults)
      .def("set_breakpenalty", &LinesP2D::set_breakpenalty)
      .def("set_error", &LinesP2D::set_error)
      .def("set_tolerance", &LinesP2D::set_tolerance)
      .def("set_verbose", &LinesP2D::set_verbose)
      .def("set_minweight", &LinesP2D::set_minweight)
      .def("set_maxoffset", &LinesP2D::set_maxoffset)
      .def("set_lsq", &LinesP2D::set_lsq)
      .def("set_unoriented", &LinesP2D::set_unoriented)
      .def("clear_ipoints", &LinesP2D::clear_ipoints)
      .def("add_ipoint", &LinesP2D::add_ipoint, py::arg("x"), py::arg("y"), py::arg("a"),
           py::arg("w") = 1.0f)
      .def("compute", py::overload_cast<>(&LinesP2D::compute))
      .def("compute", py::overload_cast<float, float, float, float>(&LinesP2D::compute))
      .def("nresults", &LinesP2D::nresults)
      .def("weight", &LinesP2D::weight)
      .def("angle", &LinesP2D::angle)
      .def("offset", &LinesP2D::offset)
      .def("nmatches", &LinesP2D::nmatches);
  m.def("makeLinesP2D", &makeLinesP2D);

  py::class_<LinesS2D>(m, "LinesS2D")
      .def("set_maxresults", &LinesS2D::set_maxresults)
      .def("set_breakpenalty", &LinesS2D::set_breakpenalty)
      .def("set_error", &LinesS2D::set_error)
      .def("set_tolerance", &LinesS2D::set_tolerance)
      .def("set_verbose", &LinesS2D::set_verbose)
      .def("set_minweight", &LinesS2D::set_minweight)
      .def("set_maxoffset", &LinesS2D::set_maxoffset)
      .def("set_lsq", &LinesS2D::set_lsq)
      .def("set_unoriented", &LinesS2D::set_unoriented)
      .def("clear_ipoints", &LinesS2D::clear_ipoints)
      .def("add_iseg", &LinesS2D::add_iseg, py::arg("x"), py::arg("y"), py::arg("x1"),
           py::arg("y1"), py::arg("a"), py::arg("w") = 1.0f)
      .def("compute", &LinesS2D::compute)
      .def("nresults", &LinesS2D::nresults)
      .def("weight", &LinesS2D::weight)
      .def("angle", &LinesS2D::angle)
      .def("offset", &LinesS2D::offset);
  m.def("makeLinesS2D", &makeLinesS2D);

  py::class_<InstanceP2D>(m, "InstanceP2D")
      .def("set_image_size", &InstanceP2D::set_image_size)
      .def("set_model_size", &InstanceP2D::set_model_size)
      .def("set_nclutter", &InstanceP2D::set_nclutter)
      .def("set_nmodel_total", &InstanceP2D::set_nmodel_total)
      .def("set_nmodel_unoccluded", &InstanceP2D::set_nmodel_unoccluded)
      .def("set_error", &InstanceP2D::set_error)
      .def("set_aerror", &InstanceP2D::set_aerror)
      .def("set_srange", &InstanceP2D::set_srange)
      .def("generate", &InstanceP2D::generate)
      .def("nimage", &InstanceP2D::nimage)
      .def("get_image",
           [](InstanceP2D &self, int i) {
             float x = 0, y = 0, a = 0;
             self.get_image(x, y, a, i);
             return py::make_tuple(x, y, a);
           })
      .def("nmodel", &InstanceP2D::nmodel)
      .def("get_model",
           [](InstanceP2D &self, int i) {
             float x = 0, y = 0, a = 0;
             self.get_model(x, y, a, i);
             return py::make_tuple(x, y, a);
           })
      .def("get_param", &InstanceP2D::get_param);
  m.def("makeInstanceP2D", &makeInstanceP2D);

  py::class_<RastP2D>(m, "RastP2D")
      .def("set_maxresults", &RastP2D::set_maxresults)
      .def("set_verbose", &RastP2D::set_verbose)
      .def("set_tolerance", &RastP2D::set_tolerance)
      .def("set_min_q", &RastP2D::set_min_q)
      .def("set_xrange", &RastP2D::set_xrange)
      .def("set_yrange", &RastP2D::set_yrange)
      .def("set_arange", &RastP2D::set_arange)
      .def("set_srange", &RastP2D::set_srange)
      .def("set_lsq", &RastP2D::set_lsq)
      .def("set_unoriented", &RastP2D::set_unoriented)
      .def("clear_msources", &RastP2D::clear_msources)
      .def("add_msource", &RastP2D::add_msource)
      .def("clear_ipoints", &RastP2D::clear_ipoints)
      .def("add_ipoint", &RastP2D::add_ipoint)
      .def("match", &RastP2D::match)
      .def("nresults", &RastP2D::nresults)
      .def("ubound", &RastP2D::ubound)
      .def("lbound", &RastP2D::lbound)
      .def("translation", &RastP2D::translation)
      .def("angle", &RastP2D::angle)
      .def("scale", &RastP2D::scale);
  m.def("makeRastP2D", &makeRastP2D);

  py::class_<RastS2D>(m, "RastS2D")
      .def("set_maxresults", &RastS2D::set_maxresults)
      .def("set_verbose", &RastS2D::set_verbose)
      .def("set_tolerance", &RastS2D::set_tolerance)
      .def("set_eps", &RastS2D::set_eps)
      .def("set_scale_eps", &RastS2D::set_scale_eps)
      .def("set_lsq", &RastS2D::set_lsq)
      .def("set_qtolerance", &RastS2D::set_qtolerance)
      .def("set_min_q", &RastS2D::set_min_q)
      .def("set_xrange", &RastS2D::set_xrange)
      .def("set_yrange", &RastS2D::set_yrange)
      .def("set_arange", &RastS2D::set_arange)
      .def("set_srange", &RastS2D::set_srange)
      .def("clear_msources", &RastS2D::clear_msources)
      .def("add_mseg", &RastS2D::add_mseg)
      .def("clear_ipoints", &RastS2D::clear_ipoints)
      .def("add_iseg", &RastS2D::add_iseg)
      .def("match", &RastS2D::match)
      .def("nresults", &RastS2D::nresults)
      .def("ubound", &RastS2D::ubound)
      .def("lbound", &RastS2D::lbound)
      .def("translation", &RastS2D::translation)
      .def("angle", &RastS2D::angle)
      .def("scale", &RastS2D::scale);
  m.def("makeRastS2D", &makeRastS2D);

  py::class_<RastSS2D>(m, "RastSS2D")
      .def("set_maxresults", &RastSS2D::set_maxresults)
      .def("set_verbose", &RastSS2D::set_verbose)
      .def("set_tolerance", &RastSS2D::set_tolerance)
      .def("set_eps", &RastSS2D::set_eps)
      .def("set_lsq", &RastSS2D::set_lsq)
      .def("set_qtolerance", &RastSS2D::set_qtolerance)
      .def("set_min_q", &RastSS2D::set_min_q)
      .def("set_xrange", &RastSS2D::set_xrange)
      .def("set_yrange", &RastSS2D::set_yrange)
      .def("set_arange", &RastSS2D::set_arange)
      .def("set_srange", &RastSS2D::set_srange)
      .def("clear_msources", &RastSS2D::clear_msources)
      .def("add_mseg", &RastSS2D::add_mseg)
      .def("clear_ipoints", &RastSS2D::clear_ipoints)
      .def("add_iseg", &RastSS2D::add_iseg)
      .def("match", &RastSS2D::match)
      .def("nresults", &RastSS2D::nresults)
      .def("ubound", &RastSS2D::ubound)
      .def("lbound", &RastSS2D::lbound)
      .def("translation", &RastSS2D::translation)
      .def("angle", &RastSS2D::angle)
      .def("scale", &RastSS2D::scale);
  m.def("makeRastSS2D", &makeRastSS2D);

  py::class_<RastRS2D>(m, "RastRS2D")
      .def("set_maxresults", &RastRS2D::set_maxresults)
      .def("set_verbose", &RastRS2D::set_verbose)
      .def("set_tolerance", &RastRS2D::set_tolerance)
      .def("set_eps", &RastRS2D::set_eps)
      .def("set_lsq", &RastRS2D::set_lsq)
      .def("set_qtolerance", &RastRS2D::set_qtolerance)
      .def("set_min_q", &RastRS2D::set_min_q)
      .def("set_xrange", &RastRS2D::set_xrange)
      .def("set_yrange", &RastRS2D::set_yrange)
      .def("set_arange", &RastRS2D::set_arange)
      .def("set_srange", &RastRS2D::set_srange)
      .def("clear_msources", &RastRS2D::clear_msources)
      .def("add_mseg", &RastRS2D::add_mseg)
      .def("clear_ipoints", &RastRS2D::clear_ipoints)
      .def("add_iseg", &RastRS2D::add_iseg)
      .def("match", &RastRS2D::match)
      .def("nresults", &RastRS2D::nresults)
      .def("ubound", &RastRS2D::ubound)
      .def("lbound", &RastRS2D::lbound)
      .def("translation", &RastRS2D::translation)
      .def("angle", &RastRS2D::angle)
      .def("scale", &RastRS2D::scale);
  m.def("makeRastRS2D", &makeRastRS2D);

  py::class_<AlignmentP2D>(m, "AlignmentP2D")
      .def("set_epsilon", &AlignmentP2D::set_epsilon)
      .def("set_srange", &AlignmentP2D::set_srange)
      .def("clear_mpoints", &AlignmentP2D::clear_mpoints)
      .def("add_mpoint", &AlignmentP2D::add_mpoint)
      .def("clear_ipoints", &AlignmentP2D::clear_ipoints)
      .def("add_ipoint", &AlignmentP2D::add_ipoint)
      .def("compute", &AlignmentP2D::compute)
      .def("quality", &AlignmentP2D::quality)
      .def("translation", &AlignmentP2D::translation)
      .def("angle", &AlignmentP2D::angle)
      .def("scale", &AlignmentP2D::scale);
  m.def("makeAlignmentP2D", &makeAlignmentP2D);

  py::class_<iupr_cedges::EdgeDetector>(m, "EdgeDetector")
      .def("set_gauss", &iupr_cedges::EdgeDetector::set_gauss)
      .def("set_noise", &iupr_cedges::EdgeDetector::set_noise)
      .def("set_poly", &iupr_cedges::EdgeDetector::set_poly)
      .def("clear", &iupr_cedges::EdgeDetector::clear)
      .def("load_pnm", &iupr_cedges::EdgeDetector::load_pnm)
      .def("save_pnm", &iupr_cedges::EdgeDetector::save_pnm)
      .def("dim", &iupr_cedges::EdgeDetector::dim)
      .def("compute", &iupr_cedges::EdgeDetector::compute)
      .def("gradient_magnitude", &iupr_cedges::EdgeDetector::gradient_magnitude)
      .def("gradient_angle", &iupr_cedges::EdgeDetector::gradient_angle)
      .def("nextchain", &iupr_cedges::EdgeDetector::nextchain)
      .def("npoints", &iupr_cedges::EdgeDetector::npoints)
      .def("point",
           [](iupr_cedges::EdgeDetector &self, int i) {
             float x = 0, y = 0;
             self.point(i, x, y);
             return py::make_tuple(x, y);
           })
      .def("nsegments", &iupr_cedges::EdgeDetector::nsegments)
      .def("segment", [](iupr_cedges::EdgeDetector &self, int i) {
        float x0 = 0, y0 = 0, x1 = 0, y1 = 0, angle = 0, magnitude = 0;
        int n = 0;
        self.segment(i, x0, y0, x1, y1, angle, magnitude, n);
        return py::make_tuple(x0, y0, x1, y1, angle, magnitude, n);
      });
  m.def("makeEdgeDetector", &iupr_cedges::makeEdgeDetector);
}
