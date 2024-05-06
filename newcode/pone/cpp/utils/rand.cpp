//
// Created by Andrei Lalu on 15/04/2024.
//

#include "rand.h"

namespace pone {
  int Rand::RandInt(int lowerBound, int upperBound)
  {
    if (upperBound > lowerBound)
      return stdUniformInt() % (upperBound - lowerBound + 1) + lowerBound;

    return (lowerBound);
  }

  double Rand::RandDouble()
  {
    return (static_cast<double>(stdUniformInt()) / RAND_MAX);
  }

  double Rand::RandDouble(double lowerBound, double upperBound)
  {
    if (upperBound > lowerBound)
      return (lowerBound + RandDouble() * (upperBound - lowerBound));

    return (lowerBound);
  }

  bool Rand::RandBool()
  {
    return (stdUniformInt() % 2 == 0);
  }

  template<typename T>
  std::vector<T> Rand::RandVector(T lowerBound, T upperBound, size_t size, bool sort)
  {
    std::vector<T> resultVec(size);
    for( size_t i=0 ; i<size ; i++)
      resultVec[i] = static_cast<T>(lowerBound + (upperBound-lowerBound) * RandDouble());

    if(sort)
      std::sort(resultVec.begin(), resultVec.end());

    return resultVec;
  }
}