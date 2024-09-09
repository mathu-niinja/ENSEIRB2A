#include <vector>
#include <string>
#include <Player.h>

int main()
{
    // Nombre de joueurs
    int nbplayers(4);

    // Les noms
    std::vector<std::string> names(nbplayers);
    // Leurs caractéristiques
    std::vector<PhysicalParameters> params(nbplayers);

    // Ines, taille moyenne, lanceuse de poids
    params[0].height = 1.63; params[0].velocity = 7; params[0].angle = 0.5;
    // Mael, taille moyenne, pas de sport connu
    params[1].height = 1.75; params[1].velocity = 6.5; params[1].angle = 0.7;
    // Abdel, taille supérieure à la moyenne, nageur
    params[2].height = 1.83; params[2].velocity = 7.3; params[2].angle = 0.6;
    // Jeanne, taille inférieure à la moyenne, pas de sport connu
    params[3].height = 1.58; params[3].velocity = 5.7; params[3].angle = 0.65;

    // Nombre de lancers
    const unsigned int nbLaunch(3);

    // Boucle sur les joueurs
    for (unsigned int p = 0; p < nbplayers; p++)
    {
        // Définir le player
        Player player(p, names[p], params[p], nbLaunch);
        // Le faire jouer
        player.throwAll();
    }

    return 0;
}