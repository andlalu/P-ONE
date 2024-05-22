//
// Created by Andrei Lalu on 15/04/2024.
//

#ifndef PONE_RAND_H
#define PONE_RAND_H

#include <cstdlib>
#include <random>

namespace pone
{
  class Rand
  {
    public:
      explicit Rand(long seed) : m_engine( seed ), m_stdUniformIntGenerator(0,RAND_MAX ), m_seed(seed ) { ; }

      //Random integer between 0 and RAND_MAX
      int RandInt() { return(stdUniformInt() ); }

      //Random integer between lower and upper bound
      int RandIntBounded(int lowerBound, int upperBound);

      //Random double between 0 and 1.0
      double RandDouble();

      //Random double between lower and upper bound
      double RandDoubleBounded(double lowerBound, double upperBound);

      //Random bool
      bool RandBool();

      // Reset random number generator
      void Reset() { m_engine.seed( m_seed ); }

      // Random vector
      template<typename T>
      std::vector<T> RandVector(T lowerBound, T upperBound, size_t size, bool sort = false);

    private:
      int stdUniformInt() { return( m_stdUniformIntGenerator(m_engine ) ); }

      std::default_random_engine m_engine;
      std::uniform_int_distribution<int> m_stdUniformIntGenerator;
      long m_seed;
  };
}

#endif //PONE_RAND_H
