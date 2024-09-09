#ifndef _FUNCTION_H

#include "DataFile.h"

class Function {
private:
   // Quelques variables privées utiles pour
   // construire la condition initiale et la solution exacte (sigma)
   const double _sigma;
   const double _xmin, _ymin;
   const double _xmax, _ymax;

public:
   // Constructeur
   Function(DataFile* data_file);

   // Condition initiale
   double InitialCondition(const double x, const double y) const;

   // Solution exacte si elle est connue (utile pour valider le code)
   double ExactSolution(const double x, const double y, const double t) const;

   // Terme source
   double SourceFunction(const double x, const double y, const double t) const;

   // Conditions limites
   double DirichletLeftBC(const double y, const double t) const;
   double DirichletRightBC(const double y, const double t) const;
   double DirichletDownBC(const double x, const double t) const;
   double DirichletUpBC(const double x, const double t) const;

   double NeumannDownBC(const double x, const double t) const;
   double NeumannUpBC(const double x, const double t) const;
   double NeumannLeftBC(const double y, const double t) const;
   double NeumannRightBC(const double y, const double t) const;
};

#define _FUNCTION_H
#endif
