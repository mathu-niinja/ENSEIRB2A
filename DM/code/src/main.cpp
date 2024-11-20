#include "TimeScheme.h"
#include <string>
#include <iostream>
#include <cmath>

#include <chrono>
#include <random>

using namespace std;
using namespace Eigen;

int main()
{
	// Définition du temps initial, du temps final et du pas de temps
	double t0(0.), tfinal(50.), dt(0.01), d(1.);
	// Définition du nombre d'itérations à partir du temps final et du pas de temps
	int nb_iterations = int(ceil(tfinal/dt));
	// Recalcul de dt
	dt = tfinal / nb_iterations;
	// Nom du fichier solution
	string results = "solution.txt";
	string str_method;

	// Définition du vecteur initial et sol exacte finale(vecteur de Eigen)
	VectorXd sol0, exactSol, omega;
	// Définition d'un pointeur de OdeSystem (nouveau pointeur vers une classe fille de OdeSystem)
	int choix_probleme, choix_schema;
	OdeSystem* sys(0);
	TimeScheme* time_scheme(0); //definition du pointeur vers la classe

    // Demander à l'utilisateur d'entrer une valeur pour l'EDO
    cout << "------------------------------------" << endl;
    cout << "Choisir l'EDO à traiter : " << endl;
    cout << "1 : 1 Vortex actif, 1 Vortex Passif" << endl;
	cout << "2 : 2 Vortex actifs avec w1 = -w2" << endl;
	cout << "3 : 2 Vortex actifs avec calcul de la solution exacte" << endl;
	cout << "4 : 4 Vortex actifs, premier exemple (choix du lambda)" << endl;
	cout << "5 : 4 Vortex actifs, deuxième exemple" << endl;
	cout << "6 : N Vortex actifs placés en ellipse" << endl;
	cout << "7 : N Vortex actifs, M Vortex passifs placés aléatoirement" << endl;
	cout << "------------------------------------" << endl;
    cin >> choix_probleme;


    // Utilisation du switch pour choisir le probleme
    switch (choix_probleme) {
		case 1: { //cas avec N=1 vortex (1 actif, 1 passif)

			int N = 2;

			omega.resize(N);
			omega(0) = 1;
			omega(1) = 0;

			sol0.resize(2*N);
			sol0(0) = 0.5;
			sol0(1) = 0;
			sol0(2) = -0.5;
			sol0(3) = 0;

			results = "N_1_Vortex";
			
			sys = new VortexSystem(N,omega);
			break; 
		}

		case 2: { //cas avec N=2 vortex, w1 = - w2

			int N = 2;
			int k = 1;
			d = 1;

			omega.resize(N);
			omega(0) = -k;
			omega(1) = k;

			sol0.resize(2*N);
			sol0(0) = d/2;
			sol0(1) = 0.0;
			sol0(2) = -d/2;
			sol0(3) = 0;

			results = "N_2_Vortex";
			
			sys = new VortexSystem(N,omega);

			/* //Bizarre la solution exacte
			// calcul de la solution exacte 
			VectorXd Cvort;
			
			exactSol.resize(2*N);
			exactSol(0) = d/2 ;
			exactSol(1) = k/(2*M_PI*d) * tfinal ;
			exactSol(2) = -d/2 ; 
			exactSol(3) = k/(2*M_PI*d) * tfinal ;
			*/

			break; 
			}
		
		case 3: { //cas avec N=2 vortex qui se tournent autour

			int N = 2;
			d = 1.;

			omega.resize(N);
			omega(0) = 1.0;
			//omega(1) = 2.0;
			omega(1) = 1.0;

			sol0.resize(2*N);
			sol0(0) = d/2;
			sol0(1) = 0.0;
			sol0(2) = -d/2;
			sol0(3) = 0;

			results = "N_2_Vortex";
			
			sys = new VortexSystem(N,omega);

			// calcul de la solution exacte 
			double gamma, theta;
			VectorXd Cvort;
			
			exactSol.resize(2*N);
			exactSol = sol0;

			Cvort.resize(2);
			Cvort(0) = (omega(0)*exactSol(0) + omega(1)*exactSol(2)) / (omega(0)+omega(1)) ; 
			Cvort(1) = (omega(0)*exactSol(1) + omega(1)*exactSol(3)) / (omega(0)+omega(1)) ; 
			gamma = omega(0) / (omega(0)+omega(1)) ;
			theta = (omega(0)+omega(1)) / (2 * M_PI * d * d) ;
			exactSol(0) = Cvort(0) + (1.-gamma)*d*cos(theta*tfinal) ;
			exactSol(1) = Cvort(1) + (1.-gamma)*d*sin(theta*tfinal) ;
			exactSol(2) = Cvort(0) - gamma * d * cos(theta*tfinal) ; 
			exactSol(3) = Cvort(1) - gamma * d * sin(theta*tfinal) ;

			break; 
			}

		case 4: { //cas avec N=4 vortex exemple 1

			int N = 4;
			int lambda;

			cout << "Choisir lambda : " << endl;
			cout << "1 : lambda = 1" << endl;
			cout << "1 : lambda = 2" << endl;
			cout << "1 : lambda = -1" << endl;

            cin >> lambda;

			omega.resize(N);
			omega(0) = 1;
			omega(1) = lambda;
			omega(2) = lambda;
			omega(3) = 1;

			sol0.resize(2*N);
			sol0(0) = 1;
			sol0(1) = 0;
			sol0(2) = -1;
			sol0(3) = 0;
			sol0(4) = -1.2;
			sol0(5) = 0;
			sol0(6) = 1.2;
			sol0(7) = 0;

			results = "N_4_Vortex";
			
			sys = new VortexSystem(N,omega);
			break; 
			}
		
		case 5: { //cas avec N=4 vortex exemple 2

			int N = 4;

			omega.resize(N);
			omega(0) = 1;
			omega(1) = 1;
			omega(2) = -1;
			omega(3) = -1;

			sol0.resize(2*N);
			sol0(0) = -1;
			sol0(1) = 0.1;
			sol0(2) = 1;
			sol0(3) = -0.1;
			sol0(4) = 1;
			sol0(5) = 0.1;
			sol0(6) = -1;
			sol0(7) = -0.1;

			results = "N_4_Vortex";
			
			sys = new VortexSystem(N,omega);
			break; 
			}

        case 6: { //cas avec N vortex

            int N;
            cout << "Choisir le nombre de vortex : " << endl;
            cin >> N;

			omega.resize(N);
			sol0.resize(2*N);
			
            for (int j = 0; j < N; j++)
            {
                omega(j) = 1;
                sol0(2*j) = cos(2*j*M_PI / N);
                sol0(2*j + 1) = sin(2*j*M_PI / N);
            }

			results = "N_" + to_string(N) + "_Vortex";
			
			sys = new VortexSystem(N,omega);
			break; 
			}

		case 7: { //cas avec N vortex actifs et M passifs

            int N;
            cout << "Choisir le nombre de vortex actifs: " << endl;
            cin >> N;

			int M;
            cout << "Choisir le nombre de vortex passifs: " << endl;
            cin >> M;

			omega.resize(N+M);
			omega.setZero();
			sol0.resize(2*(N+M));
			
			std::random_device rd;
			std::mt19937 gen(rd());
			std::uniform_real_distribution<> dist_pos(-3.0, 3.0);

			// Initialisation des vortex actifs
			for (int i = 0; i < N; i++)
			{
				sol0(2*i) = dist_pos(gen);
				sol0(2*i+1) = dist_pos(gen);
				omega(i) = 1.0; // Circulation pour vortex actif
			}

			// Initialisation des vortex passifs
			for (int i = N; i < N+M; i++)
			{
				sol0(2*i) = dist_pos(gen);
				sol0(2*i+1) = dist_pos(gen);
				omega(i) = 0.0; // Circulation pour vortex passif
			}
			

			results = "N_" + to_string(N) + "_VortexActifs_M" + to_string(M) + "_VortexPassifs";
			
			sys = new VortexSystem(N+M,omega);
			break; 
			}
        default:{
            cout << "Vous n'avez pas cette possibilité" << endl;
			return 0;}
    }
	
	// Demander à l'utilisateur d'entrer une valeur pour le schema
	cout << "------------------------------------" << endl;
    cout << "Choisir le schema à utiliser : " << endl;
	cout << "1 : Euler Explicite " << endl;
	cout << "2 : Runge Kutta 4 " << endl;
	cout << "------------------------------------" << endl;
    cin >> choix_schema;

	switch (choix_schema) {
        case 1:
		{
			time_scheme = new EulerScheme(); // Objet de Euler
			results += "_Euler.txt";
			str_method = "Euler explicite";
			break;
		}
		case 2:
		{
			time_scheme = new RungeKuttaScheme4(); // Objet de RK4
			results += "_RK4.txt";
			str_method = "Runge Kutta 4";
			break;
		}
		default:{
            cout << "Vous n'avez pas cette possibilité pour le schéma" << endl;
			return 0;}
    }
	
	auto start = chrono::high_resolution_clock::now(); //init chrono 

	time_scheme->Initialize(t0, dt, sol0, results, sys); // Initialisation
	time_scheme->SaveSolution(0); // Sauvegarde condition initiale, commenter cette ligne pour une mesure du temps plus précise 

	for (int n = 0; n < nb_iterations; n++)
	{ // Boucle en temps
		time_scheme->Advance();
		time_scheme->SaveSolution(n+1); //commenter cette ligne pour une mesure du temps plus précise
	}

	auto end = chrono::high_resolution_clock::now(); //fin chrono
	chrono::duration<double> duration = end - start;
	cout << " " << endl; 
	cout << "Temps d'exécution pour le schéma " << str_method << " = " << duration.count() << " secondes\n" <<endl;

	if ((choix_probleme == 2) || (choix_probleme == 3)){
		VectorXd approxSol = time_scheme->GetIterateSolution(); // Temps final
		double error = ((approxSol-exactSol).array().abs()).sum();
		cout << "Erreur = " << error<< " pour dt = " << dt << endl;
		
		time_scheme->Initialize(t0, dt/2., sol0, results, sys);
		for (int n = 0; n < nb_iterations*2; n++){
			time_scheme->Advance();
		}
		approxSol = time_scheme->GetIterateSolution(); // Temps final
		double error2 = ((approxSol-exactSol).array().abs()).sum();
		cout << "Erreur = " << error2 << " pour dt = " << dt/2. << endl;
		cout << "Ordre de la méthode '" << str_method << "' = " << log2(error/error2) << endl;
	}

	delete sys;
	delete time_scheme; 
	 

	return 0;
}
