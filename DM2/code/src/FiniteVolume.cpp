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
    //this->_BC_RHS.setZero();
	this->_BC_RHS.setConstant(1.);
    vector<Triplet<double>> triplets;	triplets.clear();

	string Flux_num_choice(_df->Get_numerical_flux_choice());

	double mu = this->_df->Get_mu();  // Coefficient de diffusion

	double alpha, beta, val_B_alpha, val_B_beta; 
	double d_ij,d_ie;
	double c_ab, Pe; 

    for (unsigned int k = 0; k < this->_msh->Get_edges().size(); k++)
    {
		double nx(_msh->Get_edges_normal()(k,0)),ny(_msh->Get_edges_normal()(k,1));

		Eigen::Vector2d X;

		X(0)= this->_msh->Get_edges_center()(k,0);
		X(1)= this->_msh->Get_edges_center()(k,1);
		double vn(nx*_fct->Velocity(X,t)(0)+ny*_fct->Velocity(X,t)(1));
		
		int i = _msh->Get_edges()[k].Get_T1(); //t1
		int j = _msh->Get_edges()[k].Get_T2(); //t2

		double e_k = this->_msh->Get_edges_length()[k]; //longueur arete k 
		double area_Ti = this->_msh->Get_triangles_area()[i]; //aire de t1

        if (j != (-1)) { //arrete interieur      

            X[0]=_msh->Get_triangles_center()(i,0)-_msh->Get_triangles_center()(j,0);
			X[1]=_msh->Get_triangles_center()(i,1)-_msh->Get_triangles_center()(j,1);			
			d_ij = sqrt(X(0)*X(0)+X(1)*X(1));

			// Calcul des coefficients alpha et beta pour des conditions de bord
			c_ab = mu / d_ij;
			Pe = (d_ij * vn) / mu; //vn est deja le produit de la normal et de la vitesse 

			if (Flux_num_choice=="upwind"){
				val_B_alpha = 1. + max(Pe,0.);
				val_B_beta = 1. + max(-Pe,0.);
			}
			else if (Flux_num_choice=="centered"){
				val_B_alpha = 1 + Pe/2;
				val_B_beta = 1 - Pe/2;	
			} 

			else if (Flux_num_choice=="SG"){
				if (Pe < 1e-3) {
					val_B_alpha = 1 + Pe/2;
					val_B_beta = 1 - Pe/2;
				}

				else {
					val_B_alpha = -Pe/(exp(-Pe)-1);
					val_B_beta = Pe/(exp(Pe)-1);
				}
			} 

			else {
				cout << "Choix de flux impossible" << endl; 
			}

			alpha = c_ab * val_B_alpha;
			beta = c_ab * -val_B_beta;

			// cout << "alpha peclet= " << alpha << endl;
			// cout << "BETA peclet= " << beta << endl;

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
			
			// Calcul des coefficients alpha et beta pour des conditions de bord
			c_ab = mu / d_ie;
			Pe = (d_ie * vn) / mu; //vn est deja le produit de la normal et de la vitesse 
			

			if (Flux_num_choice=="upwind"){
				val_B_alpha = 1. + max(Pe,0.);
				val_B_beta = 1. + max(-Pe,0.);
			}
			else if (Flux_num_choice=="centered"){
				val_B_alpha = 1 + Pe/2;
				val_B_beta = 1 - Pe/2;	
			} 

			else if (Flux_num_choice=="SG"){
				if (Pe < 1e-3) {
					val_B_alpha = 1 + Pe/2;
					val_B_beta = 1 - Pe/2;
				}

				else {
					val_B_alpha = -Pe/(exp(-Pe)-1);
					val_B_beta = Pe/(exp(Pe)-1);
				}
			} 

			else {
				cout << "Choix de flux impossible" << endl; 
			}

			alpha = c_ab * val_B_alpha;
			beta = c_ab * val_B_beta;

			// cout << "alpha peclet= " << alpha << endl;
			// cout << "BETA peclet= " << beta << endl;

			Eigen::Vector2d Vect_center;
			Vect_center[0]=_msh->Get_edges_center()(k,0);
			Vect_center[1]=_msh->Get_edges_center()(k,1);	

			// Gestion des conditions aux bords
			string bc_type = _msh->Get_edges()[k].Get_BC(); // Type de condition aux bords
			if (bc_type == "Dirichlet") {
				const double h = _fct->Dirichlet_value(Vect_center,t);
				triplets.push_back({i,i, alpha*e_k/area_Ti});
				_BC_RHS[i] += beta*h*e_k/area_Ti;
			}
			else if (bc_type == "Flux") {
				const double f = _fct->Flux_value(Vect_center,t);
				_BC_RHS[i] += -f*e_k/area_Ti;
			}
			else {
				std::cout << "CL " << bc_type << " inconnue!" << std::endl; exit(-1);
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
void FiniteVolume::Save_sol(const Eigen::VectorXd& sol, int n, const double t, std::string st)
{
    double integral = 0.0;
    double norm_l1 = 0.0, norm_l2 = 0.0;
    double min = std::numeric_limits<double>::max();
    double max = std::numeric_limits<double>::lowest();
    for (unsigned int i = 0; i < sol.rows(); i++)
    {
        const double A = this->_msh->Get_triangles_area()[i];
        integral += sol(i)*A;
        norm_l1 += fabs(sol(i))*A;
        norm_l2 += sol(i)*sol(i)*A;
        min = fmin(min,sol(i));
        max = fmax(max,sol(i));
    }
    norm_l2 = sqrt(norm_l2);

    if (st == "solution")
    {
        cout << "Stats for solution: " << endl;
        cout << "    integral=" << integral << " min=" << min << " max=" << max << endl;
        cout << "    l1=" << norm_l1 << " l2=" << norm_l2 << endl;
    }

    string name_file = this->_df->Get_results() + "/" + st + "_" + std::to_string(n) + ".vtk";
    unsigned int nb_vert = this->_msh->Get_vertices().size();
    assert(((long unsigned int)sol.size() == this->_msh->Get_triangles().size())
    && "The size of the solution vector is not the same than the number of _triangles !");

    cout << "Saving solution to: " << name_file << endl;

    ofstream solution;
    solution.open(name_file, ios::out);
    solution.precision(7);

    solution << "# vtk DataFile Version 3.0 " << endl;
    solution << "2D Unstructured Grid - Finite Volume solver" << endl;
    solution << "ASCII" << endl;
    solution << endl;
    solution << "DATASET UNSTRUCTURED_GRID" << endl;
    solution << "FIELD FieldData 2" << endl;
    solution << "TIME 1 1 double " << t << endl;
    solution << "CYCLE 1 1 int " << n << endl;
    solution << endl;

    solution << "POINTS " << nb_vert << " double " << endl;
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

    // Save the solution
    solution << "SCALARS sol double 1" << endl;
    solution << "LOOKUP_TABLE default" << endl;
    for (unsigned int i = 0 ; i < this->_msh->Get_triangles().size() ; ++i)
    {
        solution << sol[i] << endl;
    }
    solution << endl;

    // Save the velocity field
    // In Paraview: add a "Glyph" in order to visualize the vector field
    solution << "VECTORS velocity double" << endl;
    for (unsigned int i = 0 ; i < this->_msh->Get_triangles().size() ; ++i)
    {
        const Eigen::Vector2d &v = this->_fct->Velocity( this->_msh->Get_triangles_center().row(i), t );
        solution << v[0] << " " << v[1] << " 0.0" << endl;
    }
    solution << endl;

    // Save the Courant number
    solution << "SCALARS Courant double 1" << endl;
    solution << "LOOKUP_TABLE default" << endl;
    for (unsigned int i = 0 ; i < this->_msh->Get_triangles().size() ; ++i)
    {
        const double v = this->_fct->Velocity( this->_msh->Get_triangles_center().row(i), t ).norm();
        solution << this->_df->Get_dt() * v / this->_msh->Get_triangles_length()(i) << endl;
    }
    solution << endl;

    // Save the Péclet number
    if (this->_df->Get_mu() > std::numeric_limits<double>::min())
    {
        solution << "SCALARS Peclet double 1" << endl;
        solution << "LOOKUP_TABLE default" << endl;
        // To avoid strange behaviour (which appear only with Apple)
        // with Paraview when we have very small data (e-35 for example)
        for (unsigned int i = 0 ; i < this->_msh->Get_triangles().size() ; ++i)
        {
            const double v = this->_fct->Velocity( this->_msh->Get_triangles_center().row(i), t ).norm();
            solution << v / this->_df->Get_mu() * this->_msh->Get_triangles_length()(i) << endl;
        }
        solution << endl;
    }

    solution.close();
}

#define _FINITEVOLUME_CPP
#endif
