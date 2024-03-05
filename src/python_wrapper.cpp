#include <iostream>
#include "geometry_types.h"
#include "Curve.h"
#include "center_clustering_algs.h"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

auto trivialFilter = [](const Candidate &c){return true;};

PYBIND11_MODULE(klcluster,m){

    m.def("TRIVIALFILTER",trivialFilter);

    py::class_<Point>(m, "Point")
            .def(py::init<dimensions_t>())
            .def("__len__", &Point::dimensions)
            .def("__getitem__", &Point::get)
            .def("__setitem__", &Point::set)
            .def("__str__", &Point::str)
            .def("__iter__", [](Point &v) { return py::make_iterator(v.begin(), v.end()); }, py::keep_alive<0, 1>())
            //.def("__repr__", &Point::repr)
            .def_property_readonly("values", &Point::as_ndarray)
            ;

    py::class_<Curve>(m, "Curve")
            .def(py::init<py::array_t<distance_t>>())
            .def(py::init<py::array_t<distance_t>, std::string>())
            .def_property_readonly("dimensions", &Curve::dimensions)
            .def_property_readonly("complexity", &Curve::size)
            //.def("add",[](Curve& c, Point& p, double w=0.5){c.push_back(p,w);})
            //.def_property("name", &Curve::get_name, &Curve::set_name)
            .def_property_readonly("values", &Curve::as_ndarray)
            //.def_property_readonly("centroid", &Curve::centroid)
            .def("__getitem__", &Curve::operator[], py::return_value_policy::reference)
            .def("__len__", &Curve::size)
            .def("__str__", &Curve::str)
            .def("__iter__", [](Curve &v) { return py::make_iterator(v.begin(), v.end()); }, py::keep_alive<0, 1>())
            //.def("__repr__", &Curve::repr)
            .def(py::pickle(
                    [](const Curve &c) {
                        return py::make_tuple(c.as_ndarray(), c.get_name());
                    },
                    [](py::tuple t) {
                        const auto coords = t[0].cast<py::array_t<distance_t>>();
                        const auto name = t[1].cast<std::string>();
                        return Curve(coords, name);
                    } ))
            ;


    py::class_<Curves>(m,"Curves")
            .def(py::init<>())
            //.def_property_readonly("m", &Curves::get_m)
            .def("add", [](Curves& cc, Curve& c){cc.push_back(c);})
            //.def("simplify", &Curves::simplify)
            .def("__getitem__", [](Curves& cc, int i){return cc[i];}, py::return_value_policy::reference)
            //.def("__setitem__", &Curves::set)
            //.def("__add__", &Curves::operator+)
            .def("__len__", &Curves::size)
            //.def("__str__", &Curves::str)
            .def("__iter__", [](Curves &v) { return py::make_iterator(v.begin(), v.end()); }, py::keep_alive<0, 1>())
            //.def("__repr__", &Curves::repr)
            .def_property_readonly("values", [](Curves& cc){
                py::list l;
                for (const Curve &elem : cc) {
                    l.append(elem.as_ndarray());
                }
                return py::array_t<distance_t>(l);
            })
            //.def_property_readonly("dimensions", &Curves::dimensions)
            .def(py::pickle(
                    [](const Curves &c) {
                        py::list l;
                        for (const Curve &elem : c) {
                            l.append(py::make_tuple(elem.as_ndarray(), elem.get_name()));
                        }
                        return l;
                    },
                    [](py::list l) {
                        Curves result;
                        for (const auto &elem : l) {
                            const auto t = elem.cast<py::tuple>();
                            const auto coords = t[0].cast<py::array_t<distance_t>>();
                            const auto name = t[1].cast<std::string>();
                            result.emplace_back(coords, name);
                        }
                        return result;
                    } ))
            ;

    py::class_<CPoint>(m,"CPoint")
            .def_property_readonly("value",&CPoint::convert);

    py::class_<CInterval>(m,"CInterval")
            .def_property_readonly("start",&CInterval::getBegin)
            .def_property_readonly("end",&CInterval::getEnd)
            .def("values",&CInterval::as_ndarray);

    py::class_<Cluster>(m,"Cluster")
            .def_property_readonly("__len__",&Cluster::size)
            .def("center",&Cluster::getCenter)
            .def("__iter__",[](Cluster &c){return py::make_iterator(c.getMatching().begin(),c.getMatching().end());},py::keep_alive<0,1>())
            .def("__getitem__", &Cluster::operator[], py::return_value_policy::reference)
            .def("values",&Cluster::as_ndarray);

    py::class_<ClusteringResult>(m,"ClusteringResult")
            .def_property_readonly("__len__",&ClusteringResult::size)
            .def("__iter__",[](ClusteringResult &c){return py::make_iterator(c.begin(),c.end());},py::keep_alive<0,1>())
            .def("__getitem__", &ClusteringResult::get, py::return_value_policy::reference)
            .def("__len__",&ClusteringResult::len);

    py::class_<CurveClusterer>(m,"CurveClusterer")
            .def(py::init<>())
            .def("initCurves",&CurveClusterer::initCurves)
            .def("test",&CurveClusterer::test)
            .def("greedyCover",[](CurveClusterer& cc, int l, int rounds){return cc.greedyCover(l,rounds,trivialFilter);})
            .def("greedyCoverAgressiveFilter",[](CurveClusterer& cc, int l, int rounds){return cc.greedyCover(l,rounds,

                                                                                                              [=](const Candidate &a) {
                                                                                                                  bool withIsTrivial = false;
                                                                                                                  bool withIsDown = false;
                                                                                                                  bool istrivial = l == 1;
                                                                                                                  bool isdown = a.getEnd() < a.getBegin();
                                                                                                                  bool nontrivial_length =
                                                                                                                          cc.simplifiedCurves[a.getCurveIndex()].subcurve_length(a.getBegin(), a.getEnd()) >
                                                                                                                          2 * (1.0 + 1.0 + 2 * (7.0 / 3.0)) * cc.getDelta();
                                                                                                                  bool nontrivial_complexity = a.getEnd().getPoint() > a.getBegin().getPoint() + l / 4;
                                                                                                                  return (withIsTrivial && istrivial) || (withIsDown && isdown) || (nontrivial_length && nontrivial_complexity);});});


}