#ifndef _LAPLACIAN_H

#include <string>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "Function.h"

class Laplacian
{
private:
	Function* _fct;
	DataFile* _df;
	Eigen::VectorXd _sol, _f;
	Eigen::MatrixXd _laplacian_matrix;

public:
	// Constructeur
	Laplacian(Function* function, DataFile* data_file);
	void InitialCondition();
	void BuildLaplacianMatrix();
	void BuildSourceTerm(double t);
    

};



#define _LAPLACIAN_H
#endif
