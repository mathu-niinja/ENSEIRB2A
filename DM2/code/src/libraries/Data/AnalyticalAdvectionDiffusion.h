#ifndef _ADVECTION_DIFFUSION_H

#include "DataFile.h"
#include <Eigen/Dense>

namespace AdvectionDiffusion
{

   /**
    * Cette classe mère inclut une série de fonctions analytiques relatives
    * à l'équation d'advection-diffusion.
    * Par défaut, elles renvoient toutes "zéro".
    * Chaque classe dérivée doit réimplémenter uniquement les méthodes souhaitées.
    */
   class Analytical {
   private:
      // Pointeur de la classe DataFile pour récupérer toutes les
      // valeurs de paramètres
      const DataFile* _df;
   protected:
      // Le coefficient de diffusion
      const double _mu;

   public: // Méthodes et opérateurs de la classe
      Analytical(DataFile* data_file);
      virtual ~Analytical() {};
      virtual double Exact_solution(const Eigen::VectorXd &X, const double t) const { return 0; };
      virtual double Initial_condition(const Eigen::VectorXd &X) const { return this->Exact_solution(X,0); };
      virtual double Source_term(const Eigen::VectorXd &X, const double t) const { return 0; };
      virtual Eigen::VectorXd Velocity(const Eigen::VectorXd &X, const double t) const;
      virtual double Flux_value(const Eigen::VectorXd &X, const double t) const { return 0; };
      virtual double Dirichlet_value(const Eigen::VectorXd &X, const double t) const { return this->Exact_solution(X,t); };
   };

   class AdvectionDiffusionAll: public Analytical {
   public:
      AdvectionDiffusionAll(DataFile* data_file): Analytical(data_file) {}
      virtual double Exact_solution(const Eigen::VectorXd &X, const double t) const;
      virtual Eigen::VectorXd Velocity(const Eigen::VectorXd &X, const double t) const;
      virtual double Source_term(const Eigen::VectorXd &X, const double t) const;
      //virtual double Flux_Value(const Eigen::VectorXd &X, const double t) const;
   };

}

#define _ADVECTION_DIFFUSION_H
#endif
