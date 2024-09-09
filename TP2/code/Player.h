#include <string>
#include <vector>
#include <random>

// La classe stockant les paramètres physiques
class PhysicalParameters
{
public:
    PhysicalParameters() {}

    double height;   // La taille du joueur
    double velocity; // L'amplitude (moyenne) de la vitesse initiale
    double angle;    // L'angle (moyen) du lancer
};

// La classe définissant le joueur
class Player
{
private: // Les attributs de la classe
    // Dépendant du joueur
    // Son identifiant
    const int _id;
    // Son nom
    const std::string _name;
    // Ses paramètres physiques
    const PhysicalParameters _params;
    // Le nombre de poids lancés à chaque fois;
    const int _nbLaunch;
    // Les distances de lancer des poids
    std::vector<double> _positions;
    // Masse moyenne d'un poids
    const double _massShot;
    // Pour générer de l'aléatoire
    std::default_random_engine _generator;
    std::normal_distribution<double> _distribution;

public: // Méthodes et opérateurs de la classe
    // Constructeur qui permet de construire un objet de la classe
    Player(const int id, const std::string & name, const PhysicalParameters &params, const int nbLaunch);
    // L'action principale du joueur
    void throwAll();
    // Distance atteinte par le joueur
    const std::vector<double> &getFinalPos() const { return _positions; };
};