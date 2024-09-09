#include <iostream>
// Pour utiliser les fonctions racine (sqrt) et puissance (pow)
#include <cmath>

using namespace std;

int main()
{
    cout.precision(15);             // pour afficher une précision à 15 chiffres
    double a(2.125878159992178711); // ou double a = 2.125878159992178711;
    float f(a);                     // ou float f = a;
    cout << "En double " << a << " versus en float " << f << endl;

    string nomFonction("racine");
    cout << nomFonction << "(" << a << ") = " << sqrt(a) << endl;

    int n(3); // ou int n = 3;
    cout << a << "^" << n << " = " << pow(a, n) << endl;
    return 0;
}