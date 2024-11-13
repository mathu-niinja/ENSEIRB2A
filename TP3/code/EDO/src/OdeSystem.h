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
    void SaveSolution(const double t, const Eigen::VectorXd & sol);
    // Pour récupérer _f
    const Eigen::VectorXd &  GetF()  const;
    // Pour construire _f en fonction de votre système
    virtual void BuildF(const double t, const Eigen::VectorXd & sol) = 0;
};

// Classe fille publique d'OdeSystem : exemple 1
class FirstExampleOdeSystem : public OdeSystem //on ne met pas le mot cle virtual ici car on est dans la classe fille 
{
  public:
    void BuildF(const double t, const Eigen::VectorXd & sol); //f(X,t) = X
};

// Classe fille publique d'OdeSystem : exemple 2
class SecondExampleOdeSystem : public OdeSystem
{
  public:
    void BuildF(const double t, const Eigen::VectorXd & sol); //x' = -y , y' = x
};

// Classe fille publique d'OdeSystem : exemple 3
class ThirdExampleOdeSystem : public OdeSystem
{
  public:
    void BuildF(const double t, const Eigen::VectorXd & sol); //x' = t*x**2
};

// Classe fille publique d'OdeSystem : exemple 4 -> lokta-volterra
class LotkaVolterraOdeSystem : public OdeSystem
{
  private:
    double _a, _b, _c, _d;
  public:
    //Constructeur de la classe 
    LotkaVolterraOdeSystem(double a, double b, double c, double d);
    //construire f
    void BuildF(const double t, const Eigen::VectorXd & sol); //x' = t*x**2
};

// Classe fille publique d'OdeSystem : exemple 4 -> lokta-volterra
class PendulumOdeSystem : public OdeSystem
{
  private:
    double _l, _m, _k;
    double _g = 9.81;
  public:
    //Constructeur de la classe 
    PendulumOdeSystem(double m, double l); // _k=0
    PendulumOdeSystem(double m, double l, double k);
    //construire f
    void BuildF(const double t, const Eigen::VectorXd & sol); 
};

#define _ODE_SYSTEM_H
#endif
