#ifndef _TIME_SCHEME_H

#include "Laplacian.h"

class TimeScheme
{
protected:
   // Pointeur vers le fichier de données
   DataFile* _df;
   // Pointeur vers le laplacien
   Laplacian* _lap;

public:
   // Constructeur par défaut
   TimeScheme(DataFile* data_file, Laplacian* lap);
   // Destructeur par défaut - Si la classe ne contient pas de destructeur par défaut
   // alors le compilateur en génère un implicitement.
   virtual ~TimeScheme();
   void SaveSol(Eigen::VectorXd sol, std::string name_sol, int n);


};

class EulerScheme : public TimeScheme
{
public:
   EulerScheme(DataFile* data_file, Laplacian* lap);
};

class ImplicitEulerScheme : public TimeScheme
{
public:
   ImplicitEulerScheme(DataFile* data_file, Laplacian* lap);
};

#define _TIME_SCHEME_H
#endif
