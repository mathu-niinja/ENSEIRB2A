/* Chargement du fichier << iostream >> qui est une bibliothèque
d'affichage de messages à l'écran dans une console */
#include <iostream> 

// Choisir parmi différentes bibliothèques celles que l'on veut
using namespace std;

/* C'est ici que commence vraiment le programme. Les programmes sont
essentiellement constitués de fonctions. Chaque fonction a un rôle et peut
appeler d'autres fonctions pour effectuer certaines actions.
Tous les programmes possèdent une fonction << main >>. C'est donc la fonction
principale. Elle renvoie toujours un entier d'où la présence du << int >>
devant le main. */
int main()
{
    /* Première ligne composée de 3 éléments qui fait quelque chose de concret :
    1. cout : affichage d'un message à l'écran
    2. "Hello world!" : le message à afficher
    3. endl : retour à la ligne dans la console. */
    cout << "Hello world!" << endl;

    /* Ce type d'instruction clôt généralement les fonctions. Ici, la fonction
    main renvoie 0 pour indiquer que tout s'est bien passé (toute valeur
    différente de 0 aurait indiqué un problème). */
    return 0;
}