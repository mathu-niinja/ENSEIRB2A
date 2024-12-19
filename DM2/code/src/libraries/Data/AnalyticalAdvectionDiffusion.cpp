#ifndef _ADVECTION_DIFFUSION_CPP

#include "AnalyticalAdvectionDiffusion.h"
#include <cmath>

namespace AdvectionDiffusion
{

   Analytical::Analytical(DataFile* data_file):
   _df(data_file), _mu(data_file->Get_mu())
   {
   }

   Eigen::VectorXd Analytical::Velocity(const Eigen::VectorXd &X, const double t) const
   {
      // Fonction qui renvoie un vecteur nul de même dimension.
      // Eigen::VectorXd res(X.size());
      Eigen::VectorXd res = Eigen::VectorXd::Zero(X.size());
      return res;
   }

   Eigen::VectorXd AdvectionDiffusionAll::Velocity(const Eigen::VectorXd &X, const double t) const
   {
      const double x(X[0]), y(X[1]);
      Eigen::VectorXd res(X.size());
      //res[0] = 0.; res[1] = 0.;
      res[0] = -(y-0.5); res[1] = (x-0.5);  //peut etre un probleme de signe
      return res;
   }

   double AdvectionDiffusionAll::Source_term(const Eigen::VectorXd &X, const double t) const
   {
      const double x(X[0]), y(X[1]);
      double x0(0.9), y0(0.1), a(50);

      return -this->_mu*(pow(a, 2)*pow(-24*t*y*(-12*t*x*y + y - y0) + (40*t*pow(y, 2) + 2)*(20*t*x*pow(y, 2) + x - x0), 2)*exp(-this->_mu*t)*exp(-a*(pow(-12*t*x*y + y - y0, 2) + pow(20*t*x*pow(y, 2) + x - x0, 2))) + pow(a, 2)*pow(80*t*x*y*(20*t*x*pow(y, 2) + x - x0) + (-24*t*x + 2)*(-12*t*x*y + y - y0), 2)*exp(-this->_mu*t)*exp(-a*(pow(-12*t*x*y + y - y0, 2) + pow(20*t*x*pow(y, 2) + x - x0, 2))) - a*(288*pow(t, 2)*pow(y, 2) + (20*t*pow(y, 2) + 1)*(40*t*pow(y, 2) + 2))*exp(-this->_mu*t)*exp(-a*(pow(-12*t*x*y + y - y0, 2) + pow(20*t*x*pow(y, 2) + x - x0, 2))) - a*(3200*pow(t, 2)*pow(x, 2)*pow(y, 2) + 80*t*x*(20*t*x*pow(y, 2) + x - x0) + (-24*t*x + 2)*(-12*t*x + 1))*exp(-this->_mu*t)*exp(-a*(pow(-12*t*x*y + y - y0, 2) + pow(20*t*x*pow(y, 2) + x - x0, 2)))) - this->_mu*exp(-this->_mu*t)*exp(-a*(pow(-12*t*x*y + y - y0, 2) + pow(20*t*x*pow(y, 2) + x - x0, 2))) + 20*a*x*pow(y, 2)*(-24*t*y*(-12*t*x*y + y - y0) + (40*t*pow(y, 2) + 2)*(20*t*x*pow(y, 2) + x - x0))*exp(-this->_mu*t)*exp(-a*(pow(-12*t*x*y + y - y0, 2) + pow(20*t*x*pow(y, 2) + x - x0, 2))) - 12*a*x*y*(80*t*x*y*(20*t*x*pow(y, 2) + x - x0) + (-24*t*x + 2)*(-12*t*x*y + y - y0))*exp(-this->_mu*t)*exp(-a*(pow(-12*t*x*y + y - y0, 2) + pow(20*t*x*pow(y, 2) + x - x0, 2))) - a*(40*x*pow(y, 2)*(20*t*x*pow(y, 2) + x - x0) - 24*x*y*(-12*t*x*y + y - y0))*exp(-this->_mu*t)*exp(-a*(pow(-12*t*x*y + y - y0, 2) + pow(20*t*x*pow(y, 2) + x - x0, 2))) + 12*x*exp(-this->_mu*t)*exp(-a*(pow(-12*t*x*y + y - y0, 2) + pow(20*t*x*pow(y, 2) + x - x0, 2))) - 20*pow(y, 2)*exp(-this->_mu*t)*exp(-a*(pow(-12*t*x*y + y - y0, 2) + pow(20*t*x*pow(y, 2) + x - x0, 2)));
   }

   double AdvectionDiffusionAll::Exact_solution(const Eigen::VectorXd &X, const double t) const
   {
      const double x(X[0]), y(X[1]);
      double x0(0.9), y0(0.1), a(50);
      double Vx = -(y-0.5);
      double Vy = x-0.5;
      return exp(-a*((x-Vx*t-x0)*(x-Vx*t-x0)+(y-Vy*t-y0)*(y-Vy*t-y0)))*exp(-this->_mu*t);
   }
}

#define _ADVECTION_DIFFUSION_CPP
#endif
