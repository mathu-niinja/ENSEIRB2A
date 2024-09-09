#ifndef _DATA_FILE_CPP

#include "DataFile.h"
#include <fstream>
#include <iostream>
#include <cmath>
#include <toml/toml.hpp>

// Constructeur
DataFile::DataFile(std::string file_name)
{
   // Lecture du fichier de données
   auto config = toml::parse(file_name);

   // Other
   const auto& other = toml::find(config, "other");
   this->_results = toml::find<std::string>(other, "results");
   this->_sigma = toml::find<double>(other, "sigma");
   system(("rm -r ./" + this->_results).c_str());
   system(("mkdir -p ./" + this->_results).c_str());

   // Time
   const auto& time = toml::find(config, "time");
   this->_t0 = toml::find<double>(time, "t0");
   this->_tfinal = toml::find<double>(time, "tfinal");
   this->_dt = toml::find<double>(time, "dt");
   this->_scheme = toml::find<std::string>(time, "scheme");

   // Space
   const auto& space = toml::find(config, "space");
   this->_xmin = toml::find<double>(space, "xmin");
   this->_xmax = toml::find<double>(space, "xmax");
   this->_hx = toml::find<double>(space, "hx");
   this->_ymin = toml::find<double>(space, "ymin");
   this->_ymax = toml::find<double>(space, "ymax");
   this->_hy = toml::find<double>(space, "hy");

   // Boundary conditions
   const auto& BC = toml::find(config, "BC");
   this->_LBC = toml::find<std::string>(BC, "LeftBoundCond");
   this->_RBC = toml::find<std::string>(BC, "RightBoundCond");
   this->_DBC = toml::find<std::string>(BC, "DownBoundCond");
   this->_UBC = toml::find<std::string>(BC, "UpBoundCond");

   std::cout << "--------------------------------------------------" << std::endl;
   std::cout << "-------------- Adapt dt, hx and hy ---------------" << std::endl;
   std::cout << "-------------- xmax = xmin + (Nx+1)*hx -----------" << std::endl;
   std::cout << "-------------- ymax = ymin + (Ny+1)*hy -----------" << std::endl;
   std::cout << "-------------- tfinal = t0 + nb_it*dt ------------" << std::endl;
   // Calcul de _Nx et adaptation de _dx pour que (xmax - xmin) = (Nx+1)*hx
   this->_Nx = int(ceil((this->_xmax-this->_xmin)/this->_hx)-1);
   this->_hx = (this->_xmax-this->_xmin)/(this->_Nx+1.);
   // Calcul de _Ny et adaptation de _dy pour que (ymax - ymin) = (Ny+1)*hy
   this->_Ny = int(ceil((this->_ymax-this->_ymin)/this->_hy)-1);
   this->_hy = (this->_ymax-this->_ymin)/(this->_Ny+1.);

   // Calcul du pas de temps en fonction de la CFL
   if (this->_scheme == "ExplicitEuler")
   {
      this->_dt = 0.95*pow(this->_hy,2)*pow(this->_hx,2)/(2*this->_sigma*(pow(this->_hx,2)+pow(this->_hy,2)));
      std::cout << "The time step is fixed with the CFL condition: dt = " << this->_dt << "." << std::endl;
   }

   // Calcul du nombre d'itérations en temps
   int nb_iterations = int(ceil((this->_tfinal-this->_t0)/this->_dt));
   // Adapter le pas de temps pour avoir _tfinal = _t0 + nb_iterations*_dt
   this->_dt = (this->_tfinal-this->_t0) / nb_iterations;
   std::cout << "--------------------------------------------------" << std::endl;
}

#define _DATA_FILE_CPP
#endif
