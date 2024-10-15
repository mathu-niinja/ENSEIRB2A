#include <string>

class Student
{
  private:
    std::string _first_name, _last_name;

  protected:
    // Variables protégées donc accessibles dans les classes filles (contrairement aux variables privées)
    int _marks_number;
    double _average;

  public:
    // Constructeur
    Student();

    // Destructeur
    // Obligatoire dès qu'on a une fonction virtuelle !
    virtual ~Student();

    // Initialiser le prénom et le nom, la moyenne à 0 et le nombre de notes à 0
    void Initialize(std::string first_name, std::string last_name);

    // Si seulement le prénom est connu
    void Initialize(std::string first_name);

    // Ajouter une nouvelle note
    // Fonction virtuelle pure => définie dans les classes filles
    virtual void AddMark(double mark) = 0; //si virtuelle pure definition specifique dasn chaque classe fille 
    //si non virtuelle pure on peut creer une generale et des specifiques. 

    // Afficher le nom de l'élève et sa note
    void Print();
};

// Étudiant(e) sérieux(se)
class SeriousStudent : public Student
{
  public:
    // Fonction virtuelle pure
    void AddMark(double mark);
};

// Étudiant(e) toujours en retard
class LateStudent : public Student
{
  public:
    // Fonction virtuelle pure
    void AddMark(double mark);
};
