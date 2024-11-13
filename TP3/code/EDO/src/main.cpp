#include "TimeScheme.h"
#include <string>
#include <iostream>
#include <cmath>
using namespace std;
using namespace Eigen;

int main()
{
	// Définition du temps initial, du temps final et du pas de temps
	double t0(0.), tfinal(1.), dt(0.01);
	// Définition du nombre d'itérations à partir du temps final et du pas de temps
	int nb_iterations = int(ceil(tfinal/dt));
	// Recalcul de dt
	dt = tfinal / nb_iterations;
	// Nom du fichier solution
	string results = "solution.txt";

	// Définition du vecteur initial et sol exacte finale(vecteur de Eigen)
	VectorXd sol0, exactSol;
	// Définition d'un pointeur de OdeSystem (nouveau pointeur vers une classe fille de OdeSystem)
	int choix_probleme, choix_schema;
	OdeSystem* sys(0);
	TimeScheme* time_scheme(0); //definition du pointeur vers la classe

    // Demander à l'utilisateur d'entrer une valeur pour l'EDO
    cout << "Choisir l'EDO à traiter : ";
    cin >> choix_probleme;

	// Demander à l'utilisateur d'entrer une valeur pour le schema
    cout << "Choisir le schema à utiliser : ";
    cin >> choix_schema;

    // Utilisation du switch pour choisir le probleme
    switch (choix_probleme) {
        case 1:
			//On definit sol0
			sol0.resize(4);
			// Initialisation des valeurs
			sol0(0) = 2.; sol0(1) = 3.1; sol0(2) = -5.1; sol0(3) = 0.1;
			// Définition de la solution exacte au temps final
			exactSol = sol0*exp(tfinal); //ici la taille de ExactSol est definie implicitement par sol0 
			results = "first_ex";

			sys = new FirstExampleOdeSystem(); //def du systeme d'edo 
            break;
			
        case 2:
			sol0.resize(2); exactSol.resize(2); sol0(0) = 1; sol0(1) = -1;
			exactSol(0) = sol0(0)*cos(tfinal)-sol0(1)*sin(tfinal);
			exactSol(1) = sol0(1)*cos(tfinal)+sol0(0)*sin(tfinal);
			/*ici la taille de exactSol doit etre definit en utilisant celle de sol0 sinon jamais def
			car utilisation des composantes et non pas du vecteur en enntier*/
			results = "second_ex";

			sys = new SecondExampleOdeSystem();
			
			/*avec la méthode de Euler xn**2 + yn**2 au temps final vaut 13.6782 soit presque x0**2 + y0**2
			dt = 0.0001;
			avec un pas de temps plus petit on obtient 13.61 soit la valeur exacte voulue. 
			*/
			break;

		case 3:
			sol0.resize(1); exactSol.resize(1); sol0(0) = -2;
			exactSol(0) = 2*sol0(0)/(2-tfinal*tfinal*sol0(0));
			results = "third_ex";

			sys = new ThirdExampleOdeSystem();
			break;

		case 4:
			double a, b, c, d;
			//a = 0.8; b = 0.4; c = 0.6; d = 0.2;
			a = 1; b = 2; c = 3; d = 4;
			sol0.resize(2); sol0(0) = 0.5 * d/c; sol0(1) = 0.5 * a/b; //si on utilise les valeurs de d/c et a/b aucune population n'évolue donc cte dans le fichier de sortie 
			results = "LotkaVolterra";

			sys = new LotkaVolterraOdeSystem(a,b,c,d);
			tfinal = 100.;
			nb_iterations = int(ceil(tfinal/dt));
			dt = tfinal / nb_iterations;
			cout << tfinal << endl;
			break;
		
        default:
            cout << "Vous n'avez pas cette possibilité" << endl;
			return 0;
    }
	
	switch (choix_schema) {
        case 1:
			time_scheme = new EulerScheme(); // Objet de Euler
			results += "_Euler.txt";
			break;
		case 2:
			time_scheme = new RungeKuttaScheme4(); // Objet de RK4
			results += "_RK4.txt";
			break;
		default:
            cout << "Vous n'avez pas cette possibilité" << endl;
			return 0;
    }
	
	
	time_scheme->Initialize(t0, dt, sol0, results, sys); // Initialisation
	time_scheme->SaveSolution(); // Sauvegarde condition initiale

	for (int n = 0; n < nb_iterations; n++)
	{ // Boucle en temps
		time_scheme->Advance();
		time_scheme->SaveSolution();
	}

	if ((choix_probleme == 1) || (choix_probleme == 2) || (choix_probleme == 3))
	{
		VectorXd approxSol = time_scheme->GetIterateSolution(); // Temps final
		double error = ((approxSol-exactSol).array().abs()).sum();
		cout << "Erreur = " << error<< " pour dt = " << dt << endl;
		time_scheme->Initialize(t0, dt/2., sol0, results, sys);
		for (int n = 0; n < nb_iterations*2; n++)
			time_scheme->Advance();

		approxSol = time_scheme->GetIterateSolution(); // Temps final
		double error2 = ((approxSol-exactSol).array().abs()).sum();
		cout << "Erreur = " << error2<< " pour dt = " << dt/2. << endl;
		cout << "Ordre de la méthode = " << log2(error/error2) << endl;
	}

	//faire les constructeurs pour le pendule 
	
	/*
	// Définition d'un objet de TimeScheme
	TimeScheme time_scheme;

	cout << "------------------------------------" << endl;
	cout << "Euler Explicite" << endl;
	cout << "------------------------------------" << endl;
	cout << "Système : X' = X avec X0 = " << endl;
	cout << sol0 << endl;
	cout << "------------------------------------" << endl;

	// Initialisation
	time_scheme.Initialize(t0, dt, sol0, results, sys);
	// On sauvegarde la solution
	time_scheme.SaveSolution();

	// Boucle en temps
	for (int n = 0; n < nb_iterations; n++)
	{
		time_scheme.Advance();
		time_scheme.SaveSolution();
	}

	// On récupère la solution approchée au temps final avec dt 
	VectorXd approxSol = time_scheme.GetIterateSolution();
	
	double error = ((approxSol-exactSol).array().abs()).sum();
	cout << "Erreur = " << error << " pour dt = " << dt << endl;
	cout << "------------------------------------" << endl;


	time_scheme.Initialize(t0, dt/2., sol0, results, sys);
	for (int n = 0; n < nb_iterations*2; n++)
	{
		time_scheme.Advance();
	}

	// On récupère la solution approchée au temps final avec dt/2
	approxSol = time_scheme.GetIterateSolution();
	double error2 = ((approxSol-exactSol).array().abs()).sum();
	cout << "Erreur = " << error2<< " pour dt = " << dt/2. << endl;
	cout << "------------------------------------" << endl;

	//On affiche l'ordre de la methode 
	cout << "Ordre de la méthode = " << log2(error/error2) << endl;
	cout << "------------------------------------" << endl;
	

	//cout << "xn**2 + yn**2 = " << pow(approxSol(0),2) + pow(approxSol(1),2) << endl; //si on veut afficher decommenter

	*/
	delete sys;

	return 0;
}
