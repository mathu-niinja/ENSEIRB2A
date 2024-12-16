#ifndef _TIME_SCHEME_CPP

#include "TimeScheme.h"
#include <iostream>

using namespace Eigen;
using namespace std;

// Constructeur par défaut (ne pas oublier de mettre votre pointeur à 0 !!)
TimeScheme::TimeScheme(DataFile* data_file, FiniteVolume* adv) :
_fin_vol(adv),_df(data_file), _sol(adv->Initial_condition()), _t(_df->Get_t0())
{
}

ExplicitEulerScheme::ExplicitEulerScheme(DataFile* data_file, FiniteVolume* adv) :
TimeScheme(data_file, adv)
{
   std::cout << "Build time scheme class: ExplicitEulerScheme." << std::endl;
   std::cout << "-------------------------------------------------" << std::endl;
   _fin_vol->Build_flux_mat_and_rhs(this->_t);
}

ImplicitEulerScheme::ImplicitEulerScheme(DataFile* data_file, FiniteVolume* adv) :
TimeScheme(data_file, adv)
{
   std::cout << "Build time scheme class: ImplicitEulerScheme." << std::endl;
   std::cout << "-------------------------------------------------" << std::endl;
}

// Destructeur (car on a des fonctions virtuelles)
TimeScheme::~TimeScheme()
{
}

// Euler Explicite
void ExplicitEulerScheme::Advance()
{
   double dt(this->_df->Get_dt());
   _fin_vol->Build_flux_mat_and_rhs(this->_t);
   this->_sol+= - dt*(_fin_vol->Get_flux_matrix()*this->_sol+_fin_vol->Get_BC_RHS());
   this->_sol+=  dt*_fin_vol->Source_term(this->_t);
   this->_t+=_df->Get_dt();
}

// Euler Implicite
void ImplicitEulerScheme::Advance()
{
   SparseMatrix<double> Id(this->_sol.size(),this->_sol.size()), M;
   SparseLU<SparseMatrix<double>> solver;
   Eigen::VectorXd b;
   Id.setIdentity();
   double dt(this->_df->Get_dt());
   _fin_vol->Build_flux_mat_and_rhs(this->_t+dt);
   M=Id+dt*_fin_vol->Get_flux_matrix();
   solver.analyzePattern(M);
   solver.factorize(M);
   b=this->_sol - dt*_fin_vol->Get_BC_RHS()+dt*_fin_vol->Source_term(this->_t+dt);
   this->_sol=solver.solve(b);
   this->_t+=_df->Get_dt();
   /*std::cout<<"___________________________"<<std::endl;
   std::cout<< _fin_vol->Get_flux_matrix()<<std::endl;
   std::cout<<"___________________________"<<std::endl;
   std::cout<< _fin_vol->Get_BC_RHS()<<std::endl;
   std::cout<<"___________________________"<<std::endl;
   std::cout<< _fin_vol->Source_term(this->_t+dt)<<std::endl;*/
}

#define _TIME_SCHEME_CPP
#endif
