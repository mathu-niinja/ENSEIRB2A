#ifndef _LAPLACIAN_CPP

#include "Laplacian.h"
#include <iostream>

using namespace std;
using namespace Eigen;

// Constructeur
Laplacian::Laplacian(Function* function, DataFile* data_file) :
_fct(function), _df(data_file)
{
}

#define _LAPLACIAN_CPP
#endif
