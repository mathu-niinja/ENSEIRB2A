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

public:
	// Constructeur
	Laplacian(Function* function, DataFile* data_file);

};

#define _LAPLACIAN_H
#endif
