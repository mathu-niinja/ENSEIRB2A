#include <iostream>
#include <fstream>

using namespace std;

int main()
{
    // Ouvre un fichier
    ifstream mon_flux("mon_fichier.txt");

    // Premiere facon de lire un fichier (ligne par ligne)
    string ligne;

    getline(mon_flux, ligne);
    // On lit la première ligne : 0 1.1 1.5
    cout << "ligne 1 : " << ligne << endl;
    getline(mon_flux, ligne);
    // On lit la seconde : 1 1.1 3
    cout << "ligne 2 : " << ligne << endl; // etc ...

    // Deuxieme facon de lire : mot par mot (comme cin)
    double nombre1, nombre2;

    mon_flux >> nombre1 >> nombre2;
    // On lit les 2 premiers mots sur la ligne 3 : 2 et 1.1
    cout << "1er mot de la ligne 3 : " << nombre1 << endl;
    cout << "2eme mot de la ligne 3 : " << nombre2 << endl;

    // Troisieme facon de lire : caractere par caractere
    char a, b;

    mon_flux.get(a);
    mon_flux.get(b);
    // On lit le premier caractère sur la ligne 3 après 1.1 : un espace
    // (donc rien ne s'affiche!)
    cout << "premier caractère après premier mot de la ligne 3 : " << a << endl;
    // On lit le caractère après l'espace : 4
    cout << "caractère suivant : " << b << endl;

    mon_flux.close();
}