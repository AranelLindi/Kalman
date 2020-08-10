#ifndef KALMAN_H
#define KALMAN_H

// Standartbibliothek
#include <cstdint>  // Integer-Types
#include <iostream> // Konsole (Debugging)

// eigener Code
#include "matrix.hpp"

template <class T>
std::ostream &operator<<(std::ostream &o, const Matrix<T> &M) // to print matrix easily
{
   //o << endl;
   for (int i = 0; i < M.rows; i++)
   {
      o << "[ ";
      for (int j = 0; j < M.cols; j++)
         o << M(i, j) << " ";
      o << "]" << '\n';
   }
   o << std::endl;
   return o;
}

#endif // kalman.h