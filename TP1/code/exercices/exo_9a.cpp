#include <vector>
#include <iostream>
#include <fstream> // Pour pouvoir écrire et lire des fichiers}
#include <string>

using namespace std;

int main()
{
    int vec_size = 10;
    vector<double> const_vec(vec_size, 1.1);
    vector<double> it_vec(vec_size);

    for (int i = 0; i < vec_size; i++)
    {
        it_vec[i] = 1.5 * (i + 1);
    }

    vector<int> int_vec(vec_size);

    for (int i = 0; i < vec_size; i++)
    {
        int_vec[i] = i;

    }

    ofstream mon_flux;                   // Contruit un objet "ofstream"
    string name_file("mon_fichier.txt"); // Le nom de mon fichier
    mon_flux.open(name_file, ios::out);  // Ouvre un fichier appelé name_file

    if (mon_flux)                        // Vérifie que le fichier est bien ouvert
    {
        for (int i = 0; i < vec_size; i++) // Remplit le fichier
        {
            mon_flux << int_vec[i] << " " << const_vec[i] << " " << it_vec[i] << " " << endl;
        }
    }
    else // Renvoie un message d'erreur si ce n'est pas le cas
    {
        cout << "ERREUR: Impossible d'ouvrir le fichier." << endl;
    }

    mon_flux.close(); // Ferme le fichier
    
    return 0;
}