#include <iostream>
#include <string>

using namespace std;

int main()
{
    string repA("Euler"), repB("Runge-Kutta");

    cout << "Quel schéma en temps souhaitez-vous utiliser ? " << endl;
    cout << "A) " << repA << endl;
    cout << "B) " << repB << endl;

    char votre_reponse;

    cout << "Votre reponse (A ou B) : ";
    cin >> votre_reponse; // Récupère la réponse de l'utilisateur

    string *reponseUtilisateur = nullptr; // Un pointeur qui pointera sur la réponse

    switch (votre_reponse)
    {
    case 'A':
        reponseUtilisateur = &repA; // Déplace le pointeur sur la réponse choisie
        break;
    case 'B':
        reponseUtilisateur = &repB;
        break;
    default:
        cout << "Ce choix n'est pas valable !" << endl;
        exit(0); // Le programme s'arrête
    }

    // Utilise le pointeur pour afficher la réponse choisie
    cout << "Vous avez choisi la reponse : " << *reponseUtilisateur << endl;
    // Le pointeur n'ayant pas été défini par un "new"
    // il n'est pas nécessaire de libérer la case mémoire

    return 0;
}