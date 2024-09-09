#include <iostream>

using namespace std;

int main()
{
    int mon_entier(23);
    
    // Affichage de l'adresse de la variable en utilisant &
    // Elle sera donnée en base 16 (d'où la présence de lettres)
    cout << "L'adresse est : " << &mon_entier << endl;
    
    return 0;
}