//
// Created by Andrei Lalu on 06/05/2024.
//
#include <pybind11/pybind11.h>
#include <../utils/rand.h>

namespace pone
{
  enum PrefixType
  {
    NotSet = -1,
    Class = 0,
    ClassMethod = 1,
    Method = 2,
  };

  static std::string pyPrefix(PrefixType type)
  {
    switch (type)
    {
      case Class:
        return "CppC_";
      case ClassMethod:
        return "";
      case Method:
        return "cppF_";
      default:
      case NotSet:
        break;
    }
    return "cpp_unsup_";
  }


  std::string className(const std::string &inputName)
  {
    return pyPrefix(Class) + inputName;
  }

  std::string classMethodName(const std::string &inputName)
  {
    return pyPrefix(ClassMethod) + inputName;
  }

  std::string MethodName(const std::string &inputName)
  {
    return pyPrefix(Method) + inputName;
  }
}

namespace py = pybind11;

void test()
{
  auto testInstance = pone::Ra  lalu
      nd(25);
}

PYBIND11_MODULE(pone, instance)
{
  py::class_<pone::Rand>(instance, pone::className("Rand").data())
      .def(py::init<long>(), py::arg("seed"), "Constructor taking a seed")
      .def(pone::classMethodName("rand_int").data(), &pone::Rand::RandInt, "Random integer")
      .def(pone::classMethodName("rand_int").data(), &pone::Rand::RandIntBounded, "Random integer with bounds")
      .def(pone::classMethodName("rand_double").data(), &pone::Rand::RandDouble, "Random double")
      .def(pone::classMethodName("rand_double_bounded").data(), &pone::Rand::RandDoubleBounded, "Random double with bounds")
      .def(pone::classMethodName("rand_bool").data(), py::overload_cast<>(&pone::Rand::RandBool), "Random bool");
};