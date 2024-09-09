#include <vector>
#include <string>
#include <iostream>
#include <algorithm> // Pour la fonction 'sort'

using namespace std;

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