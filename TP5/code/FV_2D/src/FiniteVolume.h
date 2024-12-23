#ifndef _FINITEVOLUME_H

#include <string>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "libraries/Data/AnalyticalAdvectionDiffusion.h"
#include "libraries/Mesh/Mesh2D.h"

class FiniteVolume {
private:
	// Pointeur de la classe Analytical (accès à la condition initiale,
	// le terme source, la solution exacte, la vitesse, etc.)
	const AdvectionDiffusion::Analytical* _fct;
	// Pointeur de la classe DataFile pour récupérer toutes les
	// valeurs de paramètres
	const DataFile* _df;
	// Pointeur de la classe Mesh pour récupérer toutes les
	// données concernant le maillage
	const Mesh2D* _msh;

	// Matrice des flux
	Eigen::SparseMatrix<double> _mat_flux;

	// Membre de droite pour les conditions aux bords
	Eigen::VectorXd _BC_RHS;

public:
	// Constructeur
	FiniteVolume(AdvectionDiffusion::Analytical* fct, DataFile* data_file, Mesh2D* mesh);

	// Construit la matrice des flux et le membre de droite
	void Build_flux_mat_and_rhs(const double& t);

	// Renvoie la matrice qui permet d'obtenir le flux en faisant : M * sol
	const Eigen::SparseMatrix<double>& Get_flux_matrix() const {return _mat_flux;};

	// Renvoie le vecteur M*sol = rhs
	const Eigen::VectorXd& Get_BC_RHS() const {return _BC_RHS;};

	// Condition Initiale au centre des triangles
	Eigen::VectorXd Initial_condition();

	// Terme source au centre des triangles
	Eigen::VectorXd Source_term(double t);

	// Solution exacte au centre des triangles
	Eigen::VectorXd Exact_solution(double t);

	// Sauvegarde la solution
	void Save_sol(const Eigen::VectorXd& sol, int n, std::string st);
};


#define _FINITEVOLUME_H
#endif
