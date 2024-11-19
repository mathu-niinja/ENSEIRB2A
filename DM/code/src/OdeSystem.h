#ifndef _ODE_SYSTEM_H

#include <Eigen/Dense>
#include <fstream>

class OdeSystem
{
  protected:
    // Votre fonction f
    Eigen::VectorXd _f;
  private:
    // Écriture du fichier
    std::ofstream _file_out;
  public:
    // Constructeur par défaut
    OdeSystem();
    // Destructeur par défaut
    virtual ~OdeSystem();
    // Initialiser le nom du fichier solution
    void InitializeFileName(const std::string file_name);
    // Sauvegarde la solution
    void SaveSolution(const double t, const Eigen::VectorXd & sol, const int i);
    // Pour récupérer _f
    const Eigen::VectorXd &  GetF()  const;
    // Pour construire _f en fonction de votre système
    virtual void BuildF(const double t, const Eigen::VectorXd & sol) = 0;
};

class VortexSystem : public OdeSystem 
{
  private :
    int _N;
    Eigen::VectorXd _omega;
  public :
    VortexSystem(const int N, const Eigen::VectorXd & omega);
    void BuildF(const double t, const Eigen::VectorXd & sol);
};

#define _ODE_SYSTEM_H
#endif
