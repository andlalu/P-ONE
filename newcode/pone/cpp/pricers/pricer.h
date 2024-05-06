//
// Created by Andrei Lalu on 17/04/2024.
//

#ifndef PONE_PRICER_H
#define PONE_PRICER_H

#include <string>

class Pricer
{
  Pricer() :
  m_isValid(false)
  {
    ;
  }

  virtual double price() = 0;

  // private members
  bool m_isValid;
  std::string m_message;


};
#endif //PONE_PRICER_H
