#ifndef _FINITEVOLUME_CPP

#include "FiniteVolume.h"
#include <fstream>
#include <iostream>

using namespace std;
using namespace Eigen;
using namespace AdvectionDiffusion;

// Constructeur
FiniteVolume::FiniteVolume(Analytical* function, DataFile* data_file, Mesh2D* mesh) :
_fct(function), _df(data_file), _msh(mesh)
{
	std::cout << "Build finite volume class." << std::endl;
	std::cout << "-------------------------------------------------" << std::endl;
}

 // Construit la matrice des flux 
void FiniteVolume::Build_flux_mat_and_rhs(const double& t)
{
    // Matrix
    this->_mat_flux.resize(this->_msh->Get_triangles().size(), this->_msh->Get_triangles().size());
    // RHS
    this->_BC_RHS.resize(this->_msh->Get_triangles().size());
    this->_BC_RHS.setZero();
    vector<Triplet<double>> triplets;	triplets.clear();

	string Flux_num_choice(_df->Get_numerical_flux_choice());

	double mu = this->_df->Get_mu();  // Coefficient de diffusion

	double alphaD, betaD, alphaA, betaA, alpha, beta; 
	double d_ij,d_ie;

    for (unsigned int k = 0; k < this->_msh->Get_edges().size(); k++)
    {
		double nx(_msh->Get_edges_normal()(k,0)),ny(_msh->Get_edges_normal()(k,1));

		Eigen::Vector2d X;

		X(0)= this->_msh->Get_edges_center()(k,0);
		X(1)= this->_msh->Get_edges_center()(k,1);
		double vn(nx*_fct->Velocity(X,t)(0)+ny*_fct->Velocity(X,t)(1));
		
		if (Flux_num_choice=="upwind"){ //definition coefficients advection
			if (vn>=0)
			{
				alphaA=vn;
				betaA=0;
			}
			if (vn<0)
			{
				alphaA=0;
				betaA=vn;
			}
		}
		if (Flux_num_choice=="centered"){
			alphaA=vn/2;
			betaA=vn/2;
		}
		
		int i = _msh->Get_edges()[k].Get_T1(); //t1
		int j = _msh->Get_edges()[k].Get_T2(); //t2

		double e_k = this->_msh->Get_edges_length()[k]; //longueur arete k 
		double area_Ti = this->_msh->Get_triangles_area()[i]; //aire de t1

        if (j != (-1)) { //arrete interieur      

            X[0]=_msh->Get_triangles_center()(i,0)-_msh->Get_triangles_center()(j,0);
			X[1]=_msh->Get_triangles_center()(i,1)-_msh->Get_triangles_center()(j,1);			
			d_ij = sqrt(X(0)*X(0)+X(1)*X(1));

			alphaD = mu / d_ij;
			betaD = - mu / d_ij;

			alpha = alphaA + alphaD;
			beta = betaA + betaD;

			double area_Tj = this->_msh->Get_triangles_area()[j];

            // Ajouter les triplets à la matrice
            triplets.push_back(Triplet<double>(i, i, e_k * alpha/area_Ti)); //divise par Ti
			triplets.push_back(Triplet<double>(i, j, e_k * beta/area_Ti));
			triplets.push_back(Triplet<double>(j, i, -e_k * alpha/area_Tj));
            triplets.push_back(Triplet<double>(j, j, -e_k * beta/area_Tj)); //divise par Tj
        
        }
		if (j == -1) {  // Arête de bord
			// Calculer la distance entre le centre du triangle et le centre de l'arête
			d_ie = 2*sqrt(pow((_msh->Get_triangles_center()(i,0)-_msh->Get_edges_center()(k,0)),2)+pow(_msh->Get_triangles_center()(i,1)-_msh->Get_edges_center()(k,1),2));
			
			// Calcul des coefficients alphaD et betaD pour des conditions de bord
			alphaD = mu / d_ie;
			betaD = -mu / d_ie;

			Eigen::Vector2d Vect_center;
			Vect_center[0]=_msh->Get_edges_center()(k,0);
			Vect_center[1]=_msh->Get_edges_center()(k,1);	

			// Gestion des conditions aux bords
			string bc_type = _msh->Get_edges()[k].Get_BC(); // Type de condition aux bords
			if (bc_type == "Neumann") {
				// Condition de Neumann :
				double g = _fct->Neumann_value(Vect_center, t);  // Flux normal prescrit
				this->_BC_RHS(i) = (betaA * d_ie * g * e_k) / area_Ti ;
				this->_BC_RHS(i) += -mu * (e_k * g) / area_Ti;  // Contribution directe au RHS
				triplets.push_back({i, i, e_k*(alphaA+betaA)/area_Ti});
			}
			else if (bc_type == "Dirichlet") {
				// Valeur de Dirichlet prescrite (h)
				double h = this->_fct->Dirichlet_value(Vect_center, t);  // Valeur sur l'arête
				// Ajouter la contribution pour le terme de diffusion
				this->_BC_RHS(i) = 2 * (betaA * e_k * h)/area_Ti;
				this->_BC_RHS(i) += 2 * (betaD * e_k * h)/ area_Ti;  // Contribution au vecteur RHS pour T_i
				triplets.push_back({i, i, e_k*(alphaA-betaA)/area_Ti + 2*e_k*(alphaD)/area_Ti});  // Contribution à la diagonale de T_i
			}
		}

    }
    this->_mat_flux.setFromTriplets(triplets.begin(), triplets.end());
}

// --- Déjà implémenté ---
// Construit la condition initiale au centre des triangles
VectorXd FiniteVolume::Initial_condition()
{
	VectorXd sol0(this->_msh->Get_triangles().size());

	for (unsigned int i = 0; i < this->_msh->Get_triangles().size(); i++)
	{
		sol0(i) = this->_fct->Initial_condition(this->_msh->Get_triangles_center().row(i));
	}

	return sol0;
}

// Terme source au centre des triangles
VectorXd FiniteVolume::Source_term(double t)
{
	VectorXd sourceterm(this->_msh->Get_triangles().size());

	for (unsigned int i = 0; i < this->_msh->Get_triangles().size(); i++)
	{
		sourceterm(i) = this->_fct->Source_term(this->_msh->Get_triangles_center().row(i),t);
	}

	return sourceterm;
}

// Solution exacte au centre des triangles
VectorXd FiniteVolume::Exact_solution(const double t)
{
	VectorXd exactsol(this->_msh->Get_triangles().size());

	for (unsigned int i = 0; i < this->_msh->Get_triangles().size(); i++)
	{
		exactsol(i) = this->_fct->Exact_solution(this->_msh->Get_triangles_center().row(i), t);
	}
	return exactsol;
}

// Sauvegarde la solution
void FiniteVolume::Save_sol(const Eigen::VectorXd& sol, int n, std::string st)
{
	double norm = 0;
	for (unsigned int i = 0; i < sol.rows(); i++)
	{
		norm += sol(i)*sol(i)*this->_msh->Get_triangles_area()[i];
	}
	norm = sqrt(norm);

	if (st == "solution")
	{
		cout << "Norme de u = " << norm << endl;
	}

	string name_file = this->_df->Get_results() + "/" + st + "_" + std::to_string(n) + ".vtk";
	unsigned int nb_vert = this->_msh->Get_vertices().size();
	assert(((long unsigned int)sol.size() == this->_msh->Get_triangles().size())
	&& "The size of the solution vector is not the same than the number of _triangles !");

	ofstream solution;
	solution.open(name_file, ios::out);
	solution.precision(7);

	solution << "# vtk DataFile Version 3.0 " << endl;
	solution << "2D Unstructured Grid" << endl;
	solution << "ASCII" << endl;
	solution << "DATASET UNSTRUCTURED_GRID" << endl;

	solution << "POINTS " << nb_vert << " float " << endl;
	for (unsigned int i = 0 ; i < nb_vert ; ++i)
	{
		solution << ((this->_msh->Get_vertices()[i]).Get_coor())[0] << " "
		<< ((this->_msh->Get_vertices()[i]).Get_coor())[1] << " 0." << endl;
	}
	solution << endl;

	solution << "CELLS " << this->_msh->Get_triangles().size() << " "
	<< this->_msh->Get_triangles().size()*4 << endl;
	for (unsigned int i = 0 ; i < this->_msh->Get_triangles().size() ; ++i)
	{
		solution << 3 << " " << ((this->_msh->Get_triangles()[i]).Get_vertices())[0]
		<< " " << ((this->_msh->Get_triangles()[i]).Get_vertices())[1]
		<< " " << ((this->_msh->Get_triangles()[i]).Get_vertices())[2] << endl;
	}
	solution << endl;

	solution << "CELL_TYPES " << this->_msh->Get_triangles().size() << endl;
	for (unsigned int i = 0 ; i < this->_msh->Get_triangles().size() ; ++i)
	{
		solution << 5 << endl;
	}
	solution << endl;

	solution << "CELL_DATA " << this->_msh->Get_triangles().size() << endl;
	solution << "SCALARS sol float 1" << endl;
	solution << "LOOKUP_TABLE default" << endl;
	// To avoid strange behaviour (which appear only with Apple)
	// with Paraview when we have very small data (e-35 for example)
	double eps = 1.0e-10;
	for (unsigned int i = 0 ; i < this->_msh->Get_triangles().size() ; ++i)
	{
		solution << max(eps,sol[i]) << endl;
	}
	solution << endl;

	//solution << "CELL_DATA " << this->_msh->Get_triangles().size() << endl;
	solution << "SCALARS CFL float 1" << endl;
	solution << "LOOKUP_TABLE default" << endl;
	// To avoid strange behaviour (which appear only with Apple)
	// with Paraview when we have very small data (e-35 for example)
	for (unsigned int i = 0 ; i < this->_msh->Get_triangles().size() ; ++i)
	{
		solution << max(eps,this->_df->Get_dt()*fabs(sol[i])/this->_msh->Get_triangles_length()(i)) << endl;
	}
	solution << endl;

	if (this->_df->Get_mu() > 1e-10)
	{
		solution << "SCALARS Pe float 1" << endl;
		solution << "LOOKUP_TABLE default" << endl;
		// To avoid strange behaviour (which appear only with Apple)
		// with Paraview when we have very small data (e-35 for example)
		for (unsigned int i = 0 ; i < this->_msh->Get_triangles().size() ; ++i)
		{
			solution << max(eps,this->_msh->Get_triangles_length()(i)*fabs(sol[i])/this->_df->Get_mu()) << endl;
		}
		solution << endl;
	}

	solution.close();
}

#define _FINITEVOLUME_CPP
#endif
