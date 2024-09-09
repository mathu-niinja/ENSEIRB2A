#ifndef _GEOMETRY_CPP

#include<iostream>
#include "Geometry.h"

Geometry::Geometry(Mesh* mesh) : _vertices(mesh->getVertices()),
_triangles(mesh->getTriangles()), _edgesNeumann(mesh->getNeumannEdges())
{
}

// Construction des 3 fonctions de base
// (phihat0, phihat1, phihat2)
double Geometry::phihat(int hati, Eigen::Vector2d hatX)
{
  // TODO
  return 0.;
}

// Construction des 3 gradients des fonctions de base
// (gradphihat0, gradphihat1, gradphihat2)
// (indépendant du vecteur hatX)
Eigen::Vector2d Geometry::gradphihat(int hati)
{
  Eigen::Vector2d gradphi;
  // TODO
  return gradphi;
  }

// Construction de la fonction de transformation
// du triangle de référence hatK en K
Eigen::Vector2d Geometry::FK(int refTriangleK, Eigen::Vector2d hatX)
{
  Eigen::Vector2d X;
  // TODO
  return X;
}

// Construction de la fonction de transformation
// de l'arête de référence hatE en E
Eigen::Vector2d Geometry::FE(int refEdgeE, Eigen::Vector2d hatX)
{
  Eigen::Vector2d X;
  // TODO
  return X;
}

// Construction de la jacobienne de la fonction de transformation
// du triangle de référence hatK en K
// (indépendant du vecteur hatX)
Eigen::Matrix2d Geometry::JFK(int refTriangleK)
{
  Eigen::Matrix2d X;
  // TODO
  return X;
}

// Construction de la jacobienne de la fonction de transformation
// de l'arête de référence hatE en E
// (indépendant du vecteur hatX)
Eigen::Matrix2d Geometry::JFE(int refEdgeE)
{
  Eigen::Matrix2d X;
  // TODO
  return X;
}

// Construction de la valeur absolue du déterminant de la jacobienne
// de la fonction de transformation
// du triangle de référence hatK en K
// (indépendant du vecteur hatX)
double Geometry::measK(int refTriangleK)
{
  return 0.;
}

// Construction de la valeur absolue du déterminant de la jacobienne
// de la fonction de transformation
// de l'arête de référence hatE en E
// (indépendant du vecteur hatX)
double Geometry::measE(int refEdgeE)
{
  return 0.;
}

// Construction du gradient
// (gradphi_i) avec i qui est le hati-ème sommet de K
// appliqué en F(hatX)
Eigen::Vector2d Geometry::gradphi(int hati, int refTriangleK)
{
  Eigen::Vector2d X;
  // TODO
  return X;
}

// Points et poids de quadrature de la formule du milieu (intégration 2D)
void Geometry::quadraturePointsAndWeightsMidpointFormula(Eigen::VectorXd& weights, Eigen::Matrix<double, Eigen::Dynamic, 2>& points)
{
  weights.resize(3);
  // TODO
  points.resize(3,2);
  // TODO
}

// Points et poids de quadrature de la formule du milieu (intégration 1D)
void Geometry::quadraturePointsAndWeightsSimpsonFormula(Eigen::VectorXd& weights, Eigen::Matrix<double, Eigen::Dynamic, 2>& points, int numedge)
{
  weights.resize(3);
  // TODO
  points.resize(3,2);
  if (numedge == 1)
  {
    // TODO
  }
  else if (numedge == 2)
  {
    // TODO
  }
  else if (numedge == 3)
  {
    // TODO
  }
}


#define _GEOMETRY_CPP
#endif
