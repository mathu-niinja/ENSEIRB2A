#include <vector>
#include <string>
#include <iostream>
#include <algorithm> // Pour la fonction 'sort'

using namespace std;

/*
int main()
{
    vector<int> x(10); // Un tableau de 10 entiers
    for (int i = 0; i < 10; i++)
    {
        x[i] = i * (i + 1);
    }
    int vec_size = x.size(); // Retourne la taille du tableau x
    cout << "La taille du vecteur x est : " << vec_size << endl;
    x.resize(20); // Le tableau contient 20 éléments
    cout << "La nouvelle taille du vecteur x est : " << x.size() << endl;

    for (int i = 0; i < x.size(); i++)
    {
        cout << "x[" << i << "]=" << x[i] << endl;
    }
    return 0;
}
*/

int main()
{
    vector<string> S; // un tableau vide de string
    S.push_back("Pierre"); //Rajoute "Pierre" en fin de tableau.
    S.push_back("Julien"); //Rajoute "Julien" en fin de tableau.
    S.push_back("Paul"); //Rajoute "Paul" en fin de tableau.
    sort(S.begin(),S.end());
    // Les string sont triés par ordre alphabétique
    for( int i=0 ; i<S.size() ; ++i )
    {
        cout << "Bonjour " << S[i] << endl;
    }
    return 0;
}