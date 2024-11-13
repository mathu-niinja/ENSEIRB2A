#ifndef _LAPLACIAN_CPP

#include "Laplacian.h"
#include <iostream>
#include <vector>

using namespace std;
using namespace Eigen;

// Constructeur
Laplacian::Laplacian(Function* function, DataFile* data_file) :
_fct(function), _df(data_file)
{
}

void Laplacian::InitialCondition()
{
    int Nx = _df->Get_Nx();
    int Ny = _df->Get_Ny();
    _sol.resize(Nx*Ny);

    double xmin = _df->Get_xmin();
    double ymin = _df->Get_ymin();
    double hx = _df->Get_hx();
    double hy = _df->Get_hy();

    int k=0; 
    
    for (int j = 0; j < Ny; j++){  
        for (int i = 0; i < Nx; i++){ 

        const double x_i = xmin + (i+1) * hx; 
        const double y_j = ymin + (j+1) * hy; 
        
        _sol(k) = _fct->InitialCondition(x_i, y_j);
        k++; 
    }
   }
   //cout << "Condition initiale " << _sol << endl; //decommenter si on veut affiche la solution initiale
}

void Laplacian::BuildLaplacianMatrix() {
    int Nx = _df->Get_Nx(); // Nombre de points internes en x
    int Ny = _df->Get_Ny(); // Nombre de points internes en y
    double hx = _df->Get_hx();
    double hy = _df->Get_hy();

    // Coefficients pour la matrice B (direction x) et T (direction y)
    double coeff_B_center = 2 / (hx * hx) + 2 / (hy * hy); // Centre de B
    double coeff_B_x = -1 / (hx * hx); // Voisins de B dans x
    double coeff_T = -1 / (hy * hy); // Diagonale de T

    // Dimension totale de la matrice H (Nx * Ny) x (Nx * Ny)
    int matrix_size = Nx * Ny;

    // Initialisation de la matrice creuse
    Eigen::SparseMatrix<double> H(matrix_size, matrix_size);
    std::vector<Eigen::Triplet<double>> triplets;

    // Remplissage de la matrice H
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            int index = i * Ny + j; // Conversion (i, j) en indice linéaire

            // Coefficient central (i, j)
            triplets.push_back({index, index, coeff_B_center});


            // Voisins en x (i+1, j) et (i-1, j) si dans les limites
            if (j > 0) {
                triplets.push_back({index, index - 1, coeff_B_x});
            }
            if (j < Ny - 1) {
                triplets.push_back({index, index + 1, coeff_B_x});
            }

            // Voisins en y (i, j+1) et (i, j-1) si dans les limites
            if (i > 0) {
                triplets.push_back({index, index - Ny, coeff_T});
            }
            if (i < Nx - 1) {
                triplets.push_back({index, index + Ny, coeff_T});
            }
        }
    }

    // Assemblage de la matrice H à partir des triplets
    H.setFromTriplets(triplets.begin(), triplets.end());

    // Stockage de la matrice H dans la classe
    _laplacian_matrix = H;
    cout << H << endl;
}

void Laplacian::BuildSourceTerm(double t)
{
    int Nx = _df->Get_Nx();
    int Ny = _df->Get_Ny();
    _f.resize(Nx*Ny);

    double xmin = _df->Get_xmin();
    double ymin = _df->Get_ymin();
    double hx = _df->Get_hx();
    double hy = _df->Get_hy();

    int k=0; 
    
    for (int j = 0; j < Ny; j++){  
        for (int i = 0; i < Nx; i++){ 

        const double x_i = xmin + (i+1) * hx; 
        const double y_j = ymin + (j+1) * hy; 
        
        _f(k) = _fct->SourceFunction(x_i, y_j,t);
        k++; 
    }
   }
   cout << "Terme source au temps t = " << t << " : " << _f << endl; 
}

#define _LAPLACIAN_CPP
#endif
