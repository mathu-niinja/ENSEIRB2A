#include <iostream>
#include <ctime>
#include <cstdlib>

using namespace std;

int main()
{
    const int nb_essai(10);
    srand(time(0));
    const int nombre = (rand() % 100);
    int proposition(0), est_gagnant(nb_essai + 1);
    cout << "Le nombre caché est un entier compris entre 0 et 100." << endl;
    cout << "Vous avez " << nb_essai << " essais pour le trouver." << endl;

    for (int essai = 1; essai <= nb_essai; essai++)
    {
        cout << "Essai : " << essai << " : votre proposition est : ";
        cin >> proposition;
        if (proposition == nombre)
        {
            est_gagnant = essai;
            break;
        }
        else if (proposition > nombre)
        {
            cout << "Moins" << endl;
        }
        else
        {
            cout << "Plus" << endl;
        }
    }
    if (est_gagnant <= nb_essai)
    {
        if (est_gagnant == 1)
        {
            cout << "Vous êtes un devin : vous avez trouvé du premier coup !" << endl;
        }
        else
        {
            cout << "Bravo! Vous avez trouvé en " << est_gagnant << " essais." << endl;
        }
    }

    return 0;
}