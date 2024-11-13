#ifndef _TIME_SCHEME_CPP

#include "TimeScheme.h"
#include <iostream>

using namespace Eigen;
using namespace std;

// Constructeur par défaut (ne pas oublier de mettre votre pointeur à 0 !!)
TimeScheme::TimeScheme() : _sys(0)
{
}

// Destructeur
TimeScheme::~TimeScheme()
{
}

// Initialisation de vos différentes variables
void TimeScheme::Initialize(double t0, double dt, VectorXd & sol0, string results, OdeSystem* sys)
{
   this->_dt = dt;
   this->_t = t0 ;
   this->_sol0 = sol0;
   this->_sol = sol0;
   this->_sys = sys;
   if (results.size() > 0)
   {
      this->_sys->InitializeFileName(results);
   }
}

// Schéma en temps par défaut : ici Euler Explicite
// Avancer d'un pas de temps
void EulerScheme::Advance()
{
   this->_sys->BuildF(this->_t, this->_sol);
   this->_sol += this->_dt*this->_sys->GetF();
   this->_t += this->_dt;
}

void RungeKuttaScheme4::Advance()
{
   double h;
   h = this->_dt;

   Eigen::VectorXd k1,k2,k3,k4;

   this->_sys->BuildF(this->_t, this->_sol);
   k1 = this->_sys->GetF();

   this->_sys->BuildF(this->_t + h/2, this->_sol + h/2 * k1);
   k2 = this->_sys->GetF();

   this->_sys->BuildF(this->_t + h/2, this->_sol + h/2 * k2);
   k3 = this->_sys->GetF();

   this->_sys->BuildF(this->_t + h, this->_sol + h * k3);
   k4 = this->_sys->GetF();

   this->_sol += h/6 * (k1 + 2*k2 + 2*k3 + k4);
   this->_t += this->_dt;
   
}

// Enregistre la solution : fait appel à OdeSystem car la solution
// que l'on souhaite sauvegarder peut être différente de _sol SaveSolution
// le système
void TimeScheme::SaveSolution()
{
   this->_sys->SaveSolution(this->_t, this->_sol);
}

// Renvoie _sol (pratique pour calculer l'ordre de votre méthode)
const VectorXd & TimeScheme::GetIterateSolution() const
{
   return this->_sol;
}

#define _TIME_SCHEME_CPP
#endif
