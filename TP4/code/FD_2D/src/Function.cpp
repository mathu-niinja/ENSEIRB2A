#ifndef _FUNCTION_CPP

#include "Function.h"
#include <cmath>

Function::Function(DataFile* data_file)
: _sigma(data_file->Get_sigma()),
_xmin(data_file->Get_xmin()), _ymin(data_file->Get_ymin()),
_xmax(data_file->Get_xmax()), _ymax(data_file->Get_ymax())
{
}

// Condition initiale
double Function::InitialCondition(const double x, const double y) const
{
   return ExactSolution(x, y, 0);
}

// Solution exacte
double Function::ExactSolution(const double x, const double y, const double t) const
{
   return (y-this->_ymin)*(x-this->_xmin)*sin(y-this->_ymax)*sin(2*M_PI*(x-this->_xmax))*pow(t+1,2);
}

// Terme source
double Function::SourceFunction(const double x, const double y, const double t) const
{
   return -this->_sigma*(4*pow(M_PI,2)*(-this->_xmin + x)*(-this->_ymin + y)*pow(t+1,2)*sin(2*M_PI*(-this->_xmax + x))*sin(this->_ymax - y)
   - 4*M_PI*(-this->_ymin + y)*pow(t+1,2)*sin(this->_ymax - y)*cos(2*M_PI*(-this->_xmax + x)) + (-this->_xmin + x)*(-this->_ymin + y)*pow(t+1,2)*sin(2*M_PI*(-this->_xmax + x))*sin(this->_ymax - y)
   + 2*(-this->_xmin + x)*pow(t+1,2)*sin(2*M_PI*(-this->_xmax + x))*cos(this->_ymax - y)) + (-this->_xmin + x)*(-this->_ymin + y)*(-2*t - 2)*sin(2*M_PI*(-this->_xmax + x))*sin(this->_ymax - y);
}
// Conditions limites - Dirichlet
double Function::DirichletLeftBC(const double y, const double t) const
{
   return ExactSolution(this->_xmin, y, t);
}

double Function::DirichletRightBC(const double y, const double t) const
{
   return ExactSolution(this->_xmax, y, t);
}

double Function::DirichletDownBC(const double x, const double t) const
{
   return ExactSolution(x, this->_ymin, t);
}

double Function::DirichletUpBC(const double x, const double t) const
{
   return ExactSolution(x, this->_ymax, t);
}

// Conditions limites - Neumann
double Function::NeumannLeftBC(const double y, const double t) const
{
   double x = this->_xmin;
   double diffux = -2*M_PI*(-this->_xmin + x)*(-this->_ymin + y)*pow(t+1,2)*sin(this->_ymax - y)*cos(2*M_PI*(-this->_xmax + x)) - (-this->_ymin + y)*pow(t+1,2)*sin(2*M_PI*(-this->_xmax + x))*sin(this->_ymax - y);

   return diffux;
}

double Function::NeumannRightBC(const double y, const double t) const
{
   double x = this->_xmax;
   double diffux = -2*M_PI*(-this->_xmin + x)*(-this->_ymin + y)*pow(t+1,2)*sin(this->_ymax - y)*cos(2*M_PI*(-this->_xmax + x)) - (-this->_ymin + y)*pow(t+1,2)*sin(2*M_PI*(-this->_xmax + x))*sin(this->_ymax - y);

   return diffux;
}

double Function::NeumannDownBC(const double x, const double t) const
{
   double y = this->_ymin;
   double diffuy = (-this->_xmin + x)*(-this->_ymin + y)*pow(t+1,2)*sin(2*M_PI*(-this->_xmax + x))*cos(this->_ymax - y) - (-this->_xmin + x)*pow(t+1,2)*sin(2*M_PI*(-this->_xmax + x))*sin(this->_ymax - y);

   return diffuy;
}

double Function::NeumannUpBC(const double x, const double t) const
{
   double y = this->_ymax;
   double diffuy = (-this->_xmin + x)*(-this->_ymin + y)*pow(t+1,2)*sin(2*M_PI*(-this->_xmax + x))*cos(this->_ymax - y) - (-this->_xmin + x)*pow(t+1,2)*sin(2*M_PI*(-this->_xmax + x))*sin(this->_ymax - y);

   return diffuy;
}

#define _FUNCTION_CPP
#endif
