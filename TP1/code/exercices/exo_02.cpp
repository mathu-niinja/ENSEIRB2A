#include <iostream>
// Pour utiliser les fonctions racine (sqrt) et puissance (pow)
#include <cmath>

using namespace std;

int main()
{
    cout.precision(15);             // pour afficher une précision à 15 chiffres
    double a(2.125878159992178711); // ou double a = 2.125878159992178711;
    float fl(a);                     // ou float f = a;
    cout << "En double " << a << " versus en float " << fl << endl;

    string nomFonction("racine");
    cout << nomFonction << "(" << a << ") = " << sqrt(a) << endl;

    int n(3); // ou int n = 3;
    cout << a << "^" << n << " = " << pow(a, n) << endl;

    double b(2.), c(3.), d(2.1); // (b = 2)
    b++; // Ajoute la valeur 1 : cette commande a donné son nom au C++ ! (b = 3)
    b += c; // Ajouter la valeur c à b (b = 6)
    b -= d; // Retire la valeur d à b (b = 3.9)
    cout << "Après de nombreux calculs, b vaut " << b << endl;

    double e(2.), f(3.);
    double g = e/f;
    cout << "La division de " << e << " par " << f << " vaut " << g << "." << endl;
    
    int age(0);
    cout << "Quel age avez vous" << endl;
    cin >> age;
    cout << "Vous avez " << age << "ans !" << endl;

    int& laVariable(age);
    cout << "Vous avez " << laVariable << " ans ! (par référence)" << endl;

    return 0;
}