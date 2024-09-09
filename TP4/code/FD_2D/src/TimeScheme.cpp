#ifndef _TIME_SCHEME_CPP

#include "TimeScheme.h"
#include <iostream>
#include <fstream>
#include <cmath>


using namespace Eigen;
using namespace std;

// Constructeur par défaut (ne pas oublier de mettre votre pointeur à 0 !!)
TimeScheme::TimeScheme(DataFile* data_file, Laplacian* lap) :
_df(data_file), _lap(lap)
{
}

// Destructeur
TimeScheme::~TimeScheme()
{
}

EulerScheme::EulerScheme(DataFile* data_file, Laplacian* lap) :
TimeScheme(data_file,lap)
{
}

ImplicitEulerScheme::ImplicitEulerScheme(DataFile* data_file, Laplacian* lap) :
TimeScheme(data_file,lap)
{
}

#define _TIME_SCHEME_CPP
#endif
