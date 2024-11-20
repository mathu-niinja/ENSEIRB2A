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

/* // Construit la matrice des flux 
void FiniteVolume::Build_flux_mat_and_rhs(const double& t)
{
	// Matrix
	this->_mat_flux.resize(this->_msh->Get_triangles().size(),this->_msh->Get_triangles().size());
	// RHS
	this->_BC_RHS.resize(this->_msh->Get_triangles().size());
	this->_BC_RHS.setZero();
	vector<Triplet<double>> triplets;	triplets.clear();
	for (unsigned int i = 0; i < this->_msh->Get_edges().size(); i++)
	{
		// TODO
	}
	this->_mat_flux.setFromTriplets(triplets.begin(), triplets.end());
}*/

void FiniteVolume::Build_flux_mat_and_rhs(const double& t)
{
    // Dimensionner la matrice et initialiser le vecteur RHS
    this->_mat_flux.resize(this->_msh->Get_triangles().size(), this->_msh->Get_triangles().size());
    this->_BC_RHS.resize(this->_msh->Get_triangles().size());
    this->_BC_RHS.setZero();

    vector<Triplet<double>> triplets;
    triplets.clear();

    // Parcourir toutes les arêtes
    for (unsigned int i = 0; i < this->_msh->Get_edges().size(); i++)
    {
        const Edge& edge = this->_msh->Get_edges()[i];
        double length = this->_msh->Get_edges_length()(i); // Longueur de l'arête
        Vector2d normal = this->_msh->Get_edges_normal().row(i); // Normale de l'arête

        int t1 = edge.Get_T1(); // Triangle côté 1
        int t2 = edge.Get_T2(); // Triangle côté 2 (ou -1 si c'est une arête de bord)

        // Centres des triangles
        Vector2d center_t1 = this->_msh->Get_triangles_center().row(t1);
        double area_t1 = this->_msh->Get_triangles_area()(t1);

        Vector2d velocity = this->_fct->Velocity(center_t1.x(), center_t1.y(), t); // Champ de vitesse
        double advective_flux = velocity.dot(normal); // Flux advectif normal

        double diffusivity = this->_fct->Diffusivity(center_t1.x(), center_t1.y(), t); // Coefficient de diffusion
        double distance = (t2 != -1) ? 
                          (this->_msh->Get_triangles_center().row(t2) - center_t1).norm() : 
                          (center_t1 - edge.Get_Center()).norm(); // Distance entre centres

        // Gestion des flux internes ou arêtes de bord
        if (t2 != -1) // Flux interne
        {
            // Contribuer aux coefficients α et β
            double alphaD = diffusivity / distance; // Diffusion centrée
            double betaD = -alphaD;

            double alphaA = (advective_flux >= 0) ? advective_flux : 0; // Upwind advectif
            double betaA = (advective_flux < 0) ? advective_flux : 0;

            // Ajouter contributions à A
            triplets.emplace_back(t1, t1, (length / area_t1) * (alphaD + alphaA));
            triplets.emplace_back(t1, t2, (length / area_t1) * (betaD + betaA));
        }
        else // Condition de bord
        {
            std::string BC_type = edge.Get_BC(); // Type de condition aux bords
            if (BC_type == "Neumann")
            {
                // Neumann : ajuster RHS pour diffusion ou advection
                double g = this->_fct->Neumann_BC(edge.Get_Center().x(), edge.Get_Center().y(), t);
                this->_BC_RHS(t1) += -(length / area_t1) * diffusivity * g;
            }
            else if (BC_type == "Dirichlet")
            {
                // Dirichlet : ajuster A et RHS
                double h = this->_fct->Dirichlet_BC(edge.Get_Center().x(), edge.Get_Center().y(), t);
                double alphaD = diffusivity / distance;
                double betaD = -alphaD;

                triplets.emplace_back(t1, t1, (length / area_t1) * (alphaD));
                this->_BC_RHS(t1) += (length / area_t1) * (-betaD) * h;
            }
        }
    }

    // Assembler la matrice des flux
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
