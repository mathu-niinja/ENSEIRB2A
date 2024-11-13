#ifndef _ODE_SYSTEM_CPP

#include "OdeSystem.h"
#include <iostream>

using namespace Eigen;
using namespace std;

// Constructeur par défaut
OdeSystem::OdeSystem()
{
}

// Destructeur par défaut
OdeSystem::~OdeSystem()
{
}

// Initialisation du nom du fichier
void OdeSystem::InitializeFileName(const std::string file_name)
{
  _file_out.open(file_name);
}

// Renvoie le vecteur _f
const VectorXd & OdeSystem::GetF() const
{
  return this->_f;
}

// Enregistre la solution
// Pour le moment : sol_1, sol_2 ...
void OdeSystem::SaveSolution(const double t, const VectorXd & sol)
{
  this->_file_out << t;
  for (int i = 0 ; i < sol.rows() ; ++i)
  {
    this->_file_out << " " << sol(i);
  }
  this->_file_out << std::endl;
}

//Construction de f pour le premier exemple
void FirstExampleOdeSystem::BuildF(const double t, const VectorXd & sol)
{
  _f = sol; // f(t,X) = X pour résoudre X' = X
}

//Construction de f pour le deuxieme exemple
void SecondExampleOdeSystem::BuildF(const double t, const VectorXd & sol)
{
  _f = 0. * sol;
  _f(0) = -1 * sol(1);
  _f(1) = sol(0);
}

//Construction de f pour le troisieme exemple
void ThirdExampleOdeSystem::BuildF(const double t, const VectorXd & sol)
{
  _f = 0. * sol;
  _f(0) = t * pow(sol(0),2);
}

//Constructeur classe lotka 
LotkaVolterraOdeSystem::LotkaVolterraOdeSystem(double a, double b, double c, double d)
{
  this -> _a = a;
  this -> _b = b;
  this -> _c = c;
  this -> _d = d;
}

//Construction de f pour Lotka volterra
void LotkaVolterraOdeSystem::BuildF(const double t, const VectorXd & sol)
{
  _f = 0. * sol;
  _f(0) = sol(0)*(_a - _b * sol(1));
  _f(1) = sol(1)*(_c * sol(0) - _d);
}

//Constructeurs classe pendule
PendulumOdeSystem::PendulumOdeSystem(double m, double l)
{
  this -> _m = m;
  this -> _l = l;
  this -> _k = 0;
}

PendulumOdeSystem::PendulumOdeSystem(double m, double l, double k)
{
  this -> _m = m;
  this -> _l = l;
  this -> _k = k;
}

void PendulumOdeSystem::BuildF(const double t, const VectorXd & sol)
{
  _f = 0. * sol;
  _f(0) = sol(1)*(-_k/(_m*pow(_l,2))) - sol(0)* _g/_l;
  _f(1) = sol(1);
}



#define _ODE_SYSTEM_CPP
#endif
