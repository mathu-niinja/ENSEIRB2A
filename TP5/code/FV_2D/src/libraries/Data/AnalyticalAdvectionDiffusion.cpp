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

   double AdvectionHomogeneousNeumann::Initial_condition(const Eigen::VectorXd &X) const
   {
      const double x(X[0]), y(X[1]);
      double x0(0.5), y0(0.5), a(100);
      return 0.5*exp(-a*(pow(x-x0, 2) + pow(y-y0, 2)));
   }

   Eigen::VectorXd AdvectionHomogeneousNeumann::Velocity(const Eigen::VectorXd &X, const double t) const
   {
      Eigen::VectorXd res(X.size());
      res[0] = 0.5; res[1] = -0.8;
      return res;
   }
   Eigen::VectorXd AdvectionDiffusionAll::Velocity(const Eigen::VectorXd &X, const double t) const
   {
      const double x(X[0]), y(X[1]);
      Eigen::VectorXd res(X.size());
      res[0] = -20.0*y*y*x; res[1] = 12.0*x*y;
      return res;
   }


   double DiffusionHomogeneousNeumann::Source_term(const Eigen::VectorXd &X, const double t) const
   {
      const double x(X[0]), y(X[1]);
      double xmin(0), ymin(0), xmax(1), ymax(1);
      double x0(0.5), y0(0.5), a(5);
      return 12.5*M_PI*pow(x-xmax, 2)*pow(x-xmin, 2)*pow(y-ymax, 2)*pow(y-ymin, 2)*exp(-a*(pow(-x0+x, 2)+pow(-y0+y, 2)))*cos(0.25*M_PI*t)-this->_mu*(50.0*pow(a, 2)*pow(-2*x0+2*x, 2)*pow(x-xmax, 2)*pow(x-xmin, 2)*pow(y-ymax, 2)*pow(y-ymin, 2)*(sin(0.25*M_PI*t)+1)*exp(-a*(pow(-x0+x, 2)+pow(-y0+y, 2)))+50.0*pow(a, 2)*pow(-2*y0+2*y, 2)*pow(x-xmax, 2)*pow(x-xmin, 2)*pow(y-ymax, 2)*pow(y-ymin, 2)*(sin(0.25*M_PI*t)+1)*exp(-a*(pow(-x0+x, 2)+pow(-y0+y, 2)))-100.0*a*(-2*x0+2*x)*pow(x-xmax, 2)*(2*x-2*xmin)*pow(y-ymax, 2)*pow(y-ymin, 2)*(sin(0.25*M_PI*t)+1)*exp(-a*(pow(-x0+x, 2)+pow(-y0+y, 2)))-50.0*a*(-2*x0+2*x)*pow(x-xmin, 2)*(2*x-2*xmax)*pow(y-ymax, 2)*pow(y-ymin, 2)*(sin(0.25*M_PI*t)+1)*exp(-a*(pow(-x0+x, 2)+pow(-y0+y, 2)))-a*(-2*x0+2*x)*pow(x-xmin, 2)*(100.0*x-100.0*xmax)*pow(y-ymax, 2)*pow(y-ymin, 2)*(sin(0.25*M_PI*t)+1)*exp(-a*(pow(-x0+x, 2)+pow(-y0+y, 2)))-100.0*a*(-2*y0+2*y)*pow(x-xmax, 2)*pow(x-xmin, 2)*pow(y-ymax, 2)*(2*y-2*ymin)*(sin(0.25*M_PI*t)+1)*exp(-a*(pow(-x0+x, 2)+pow(-y0+y, 2)))-100.0*a*(-2*y0+2*y)*pow(x-xmax, 2)*pow(x-xmin, 2)*pow(y-ymin, 2)*(2*y-2*ymax)*(sin(0.25*M_PI*t)+1)*exp(-a*(pow(-x0+x, 2)+pow(-y0+y, 2)))-200.0*a*pow(x-xmax, 2)*pow(x-xmin, 2)*pow(y-ymax, 2)*pow(y-ymin, 2)*(sin(0.25*M_PI*t)+1)*exp(-a*(pow(-x0+x, 2)+pow(-y0+y, 2)))+100.0*pow(x-xmax, 2)*pow(x-xmin, 2)*pow(y-ymax, 2)*(sin(0.25*M_PI*t)+1)*exp(-a*(pow(-x0+x, 2)+pow(-y0+y, 2)))+100.0*pow(x-xmax, 2)*pow(x-xmin, 2)*pow(y-ymin, 2)*(sin(0.25*M_PI*t)+1)*exp(-a*(pow(-x0+x, 2)+pow(-y0+y, 2)))+100.0*pow(x-xmax, 2)*pow(x-xmin, 2)*(2*y-2*ymax)*(2*y-2*ymin)*(sin(0.25*M_PI*t)+1)*exp(-a*(pow(-x0+x, 2)+pow(-y0+y, 2)))+100.0*pow(x-xmax, 2)*pow(y-ymax, 2)*pow(y-ymin, 2)*(sin(0.25*M_PI*t)+1)*exp(-a*(pow(-x0+x, 2)+pow(-y0+y, 2)))+100.0*pow(x-xmin, 2)*pow(y-ymax, 2)*pow(y-ymin, 2)*(sin(0.25*M_PI*t)+1)*exp(-a*(pow(-x0+x, 2)+pow(-y0+y, 2)))+2*(2*x-2*xmin)*(100.0*x-100.0*xmax)*pow(y-ymax, 2)*pow(y-ymin, 2)*(sin(0.25*M_PI*t)+1)*exp(-a*(pow(-x0+x, 2)+pow(-y0+y, 2))));
   }

   double DiffusionAll::Source_term(const Eigen::VectorXd &X, const double t) const
   {
      const double x(X[0]), y(X[1]);
      double x0(0.5), y0(0.5), a(5);
      return 0.125*M_PI*exp(-a*(pow(x-x0, 2) + pow(y-y0, 2)))*cos(0.25*M_PI*t)-this->_mu*(pow(a, 2)*pow(2*x-2*x0, 2)*(0.5*sin(0.25*M_PI*t) + 0.5)*exp(-a*(pow(x-x0, 2) + pow(y-y0, 2))) + pow(a, 2)*pow(2*y-2*y0, 2)*(0.5*sin(0.25*M_PI*t) + 0.5)*exp(-a*(pow(x-x0, 2) + pow(y-y0, 2)))-4*a*(0.5*sin(0.25*M_PI*t) + 0.5)*exp(-a*(pow(x-x0, 2) + pow(y-y0, 2))));
   }

   double AdvectionDiffusionAll::Source_term(const Eigen::VectorXd &X, const double t) const
   {
      const double x(X[0]), y(X[1]);
      double x0(0.9), y0(0.1), a(50);

      return -this->_mu*(pow(a, 2)*pow(-24*t*y*(-12*t*x*y + y - y0) + (40*t*pow(y, 2) + 2)*(20*t*x*pow(y, 2) + x - x0), 2)*exp(-this->_mu*t)*exp(-a*(pow(-12*t*x*y + y - y0, 2) + pow(20*t*x*pow(y, 2) + x - x0, 2))) + pow(a, 2)*pow(80*t*x*y*(20*t*x*pow(y, 2) + x - x0) + (-24*t*x + 2)*(-12*t*x*y + y - y0), 2)*exp(-this->_mu*t)*exp(-a*(pow(-12*t*x*y + y - y0, 2) + pow(20*t*x*pow(y, 2) + x - x0, 2))) - a*(288*pow(t, 2)*pow(y, 2) + (20*t*pow(y, 2) + 1)*(40*t*pow(y, 2) + 2))*exp(-this->_mu*t)*exp(-a*(pow(-12*t*x*y + y - y0, 2) + pow(20*t*x*pow(y, 2) + x - x0, 2))) - a*(3200*pow(t, 2)*pow(x, 2)*pow(y, 2) + 80*t*x*(20*t*x*pow(y, 2) + x - x0) + (-24*t*x + 2)*(-12*t*x + 1))*exp(-this->_mu*t)*exp(-a*(pow(-12*t*x*y + y - y0, 2) + pow(20*t*x*pow(y, 2) + x - x0, 2)))) - this->_mu*exp(-this->_mu*t)*exp(-a*(pow(-12*t*x*y + y - y0, 2) + pow(20*t*x*pow(y, 2) + x - x0, 2))) + 20*a*x*pow(y, 2)*(-24*t*y*(-12*t*x*y + y - y0) + (40*t*pow(y, 2) + 2)*(20*t*x*pow(y, 2) + x - x0))*exp(-this->_mu*t)*exp(-a*(pow(-12*t*x*y + y - y0, 2) + pow(20*t*x*pow(y, 2) + x - x0, 2))) - 12*a*x*y*(80*t*x*y*(20*t*x*pow(y, 2) + x - x0) + (-24*t*x + 2)*(-12*t*x*y + y - y0))*exp(-this->_mu*t)*exp(-a*(pow(-12*t*x*y + y - y0, 2) + pow(20*t*x*pow(y, 2) + x - x0, 2))) - a*(40*x*pow(y, 2)*(20*t*x*pow(y, 2) + x - x0) - 24*x*y*(-12*t*x*y + y - y0))*exp(-this->_mu*t)*exp(-a*(pow(-12*t*x*y + y - y0, 2) + pow(20*t*x*pow(y, 2) + x - x0, 2))) + 12*x*exp(-this->_mu*t)*exp(-a*(pow(-12*t*x*y + y - y0, 2) + pow(20*t*x*pow(y, 2) + x - x0, 2))) - 20*pow(y, 2)*exp(-this->_mu*t)*exp(-a*(pow(-12*t*x*y + y - y0, 2) + pow(20*t*x*pow(y, 2) + x - x0, 2)));
   }

   double DiffusionAll::Neumann_value(const Eigen::VectorXd &X, const double t) const
   {
      const double x(X[0]), y(X[1]);
      double x0(0.5), y0(0.5), a(5);
      double xmin(0), ymin(0), xmax(1), ymax(1);
      double diffx = -a*(-24*t*y*(-12*t*x*y + y - y0) + (40*t*pow(y, 2) + 2)*(20*t*x*pow(y, 2) + x - x0))*exp(-this->_mu*t)*exp(-a*(pow(-12*t*x*y + y - y0, 2) + pow(20*t*x*pow(y, 2) + x - x0, 2)));
      double diffy = -a*(80*t*x*y*(20*t*x*pow(y, 2) + x - x0) + (-24*t*x + 2)*(-12*t*x*y + y - y0))*exp(-this->_mu*t)*exp(-a*(pow(-12*t*x*y + y - y0, 2) + pow(20*t*x*pow(y, 2) + x - x0, 2)));
      double nx = (fabs(x-xmin) < 1e-6) ? -1 : (fabs(x-xmax) < 1e-6) ? 1 : 0;
      double ny = (fabs(y-ymin) < 1e-6) ? -1 : (fabs(y-ymax) < 1e-6) ? 1 : 0;
      return diffx*nx+diffy*ny;
   }
   
   double AdvectionDiffusionAll::Neumann_value(const Eigen::VectorXd &X, const double t) const
   {
      const double x(X[0]), y(X[1]);
      double x0(0.9), y0(0.1), a(50);
      double xmin(0), ymin(0), xmax(1), ymax(1);
      double diffx = -a*(-12*t*y*(-6*t*x*y + y - y0) + (20*t*pow(y, 2) + 2)*(10*t*x*pow(y, 2) + x - x0))*exp(-this->_mu*t)*exp(-a*(pow(-6*t*x*y + y - y0, 2) + pow(10*t*x*pow(y, 2) + x - x0, 2)));
      double diffy = -a*(40*t*x*y*(10*t*x*pow(y, 2) + x - x0) + (-12*t*x + 2)*(-6*t*x*y + y - y0))*exp(-this->_mu*t)*exp(-a*(pow(-6*t*x*y + y - y0, 2) + pow(10*t*x*pow(y, 2) + x - x0, 2)));
      double nx = (fabs(x-xmin) < 1e-6) ? -1 : (fabs(x-xmax) < 1e-6) ? 1 : 0;
      double ny = (fabs(y-ymin) < 1e-6) ? -1 : (fabs(y-ymax) < 1e-6) ? 1 : 0;
      return diffx*nx+diffy*ny;
   }

   double DiffusionHomogeneousNeumann::Exact_solution(const Eigen::VectorXd &X, const double t) const
   {
      const double x(X[0]), y(X[1]);
      double xmin(0), ymin(0), xmax(1), ymax(1);
      double x0(0.5), y0(0.5), a(5);
      return 50.0*pow(x-xmax, 2)*pow(x-xmin, 2)*pow(y-ymax, 2)*pow(y-ymin, 2)*(sin(0.25*M_PI*t)+1)*exp(-a*(pow(-x0+x, 2)+pow(-y0+y, 2)));
   }

   double DiffusionAll::Exact_solution(const Eigen::VectorXd &X, const double t) const
   {
      const double x(X[0]), y(X[1]);
      double x0(0.5), y0(0.5), a(5);
      return 0.5*(sin(0.25*M_PI*t) + 1)*exp(-a*(pow(x-x0, 2) + pow(y-y0, 2)));
   }
   
   double AdvectionHomogeneousNeumann::Exact_solution(const Eigen::VectorXd &X, const double t) const
   {
      // Caractéristique rétrograde
      return Initial_condition(X - this->Velocity(X,t)*t);
   }   

   double AdvectionDiffusionAll::Exact_solution(const Eigen::VectorXd &X, const double t) const
   {
      const double x(X[0]), y(X[1]);
      double x0(0.9), y0(0.1), a(50);
      double Vx = -20*y*y*x;
      double Vy = 12*x*y;
      return exp(-a*((x-Vx*t-x0)*(x-Vx*t-x0)+(y-Vy*t-y0)*(y-Vy*t-y0)))*exp(-this->_mu*t);
   }
   
}

#define _ADVECTION_DIFFUSION_CPP
#endif
