// Les includes dont vous avez besoin dans le fichier .h
// Attention à ne pas mettre de using namespace dans le .h car c'est ce
// fichier qui va être inclu et il peut donc contaminer le reste du code ...
#include <iostream>
#include <vector>

// Déclaration de votre classe
class Shot
{
private: // Les attributs de la classe
    // La masse du poids
    const double _mass;
    // La position initiale du poids
    const double _y0;
    // L'amplitude moyenne de la vitesse initiale
    const double _normv0;
    // L'angle moyen de la vitesse initiale
    const double _anglev0;
    // L'angle moyen de la vitesse initiale
    std::vector<double> _v0;
    // Temps, position x, position y
    std::vector<std::vector<double>> _t_x_y;

public: // Méthodes et opérateurs de la classe
    // Constructeur qui permet de construire un objet de la classe
    // Il peut y avoir différents constructeurs
    // Il doit porter le même nom que la classe !
    Shot(const double mass, const double y0, const double normv0, const double anglev0);
    // Un deuxième constructeur (dit par copie)
    Shot(const Shot &shot);
    // Lancer un poids
    void computeTrajectory();
    // Sauvegarder la solution
    void saveSolution(std::string nameSol);

    // Récupérer la masse du poids à l'extérieur de la classe
    // (sans pouvoir la modifier !!)
    const double &getMass() const { return _mass; };
    // Récupérer la position y0
    const double &gety0() const { return _y0; };
    // Récupérer la vitesse v0
    const double &getnormv0() const { return _normv0; };
    const double &getanglev0() const { return _anglev0; };
    // Récupérer le temps, la position et la vitesse du poids
    const std::vector<std::vector<double>> &getdata() const { return _t_x_y; };
};
// Permet d'afficher les propriétés du poids
std::ostream &operator<<(std::ostream &out, Shot const &shot);