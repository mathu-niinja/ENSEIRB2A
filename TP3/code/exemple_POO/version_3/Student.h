#include <string>

class Student
{
  private:
    std::string _first_name, _last_name;
    int _marks_number;
    double _average;

  public:
    // Constructeur
    Student();

    // Initialiser le prénom et le nom, la moyenne à 0 et le nombre de notes à 0
    void Initialize(std::string first_name, std::string last_name);

    // Si seulement le prénom est connu
    void Initialize(std::string first_name);

    // Ajouter une nouvelle note
    void AddMark(double mark);

    // Afficher le nom de l'élève et sa note
    void Print();
};
