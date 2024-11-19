#ifndef _ODE_SYSTEM_CPP

#include "OdeSystem.h"
#include <iostream>
#include <cmath>
#include <iomanip>
#include <sstream>


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
  //_file_out.open(file_name);
}

// Renvoie le vecteur _f
const VectorXd & OdeSystem::GetF() const
{
  return this->_f;
}
/*
// Enregistre la solution
// Pour le moment : sol_1, sol_2 ...

void OdeSystem::SaveSolution(const double t, const VectorXd & sol, const int i)
{
  this->_file_out << t;
  for (int i = 0 ; i < sol.rows() ; ++i)
  {
    this->_file_out << " " << sol(i);
  }
  this->_file_out << std::endl;
}
*/
void OdeSystem::SaveSolution(const double t, const VectorXd & sol, const int i)
{
    // Chemin du dossier de sortie
    const std::string output_dir = "results";

    // Créer le nom de fichier avec le compteur i
    std::ostringstream filename;
    filename << output_dir << "/points_vortex.txt_" << i;

    // Ouvrir un fichier en écriture avec le nom correspondant au compteur i
    std::ofstream file_out(filename.str());
    if (!file_out.is_open()) {
        std::cerr << "Erreur : Impossible d'ouvrir le fichier " << filename.str() << std::endl;
        return;
    }

    // Écrire le temps et les valeurs de la solution dans le fichier
    //file_out << t;
    int N = sol.size()/2;
    for (int j = 0; j < N ; ++j) {
        file_out << sol(2*j) << " " << sol(2*j + 1) << std::endl; 
    }
    file_out << std::endl;

    // Le fichier sera automatiquement fermé en sortant de la portée
}

VortexSystem::VortexSystem(const int N, const VectorXd & omega)
{
  this->_N = N;
  this->_omega = omega;
}

void VortexSystem::BuildF(const double t, const VectorXd & sol)
{
  Vector2d  v_rot, v_diff, v_f;
  double norm_diff;

  _f.resize(2 * _N);

  for (int i = 0; i < _N; i++)
  {
    v_f(0) = 0;
    v_f(1) = 0;
    for (int j = 0; j < _N; j++)
    {
      if (j != i) {

        v_diff(0) = sol(2*i) - sol(2*j);
        v_diff(1) = sol(2*i + 1) - sol(2*j + 1); 
        norm_diff = v_diff.norm();

        v_rot(0) = -v_diff(1);
        v_rot(1) = v_diff(0);

        v_f += _omega(j) * v_rot/pow(norm_diff,2);
      }
    }

    v_f /= 2*M_PI;
    _f(2*i) = v_f(0);
    _f(2*i + 1) = v_f(1);
  } 
}



#define _ODE_SYSTEM_CPP
#endif
