#include "Student.h"
#include <iostream>

// -------------- Étudiant(e) -------------------------------
// Constructeur par défaut
Student::Student() {}

// Destructeur par défaut
Student::~Student() {}

// Initialiser avec prénom + nom de famille
void Student::Initialize(std::string first_name, std::string last_name)
{
  this->_first_name = first_name;
  this->_last_name = last_name;
  this->_marks_number = 0;
  this->_average = 0.;
}

// Initialiser avec prénom
void Student::Initialize(std::string first_name)
{
  this->_first_name = first_name;
  this->_last_name = "";
  this->_marks_number = 0;
  this->_average = 0.;
}


// Afficher le nom de l'étudiant(e) et sa note
void Student::Print()
{
  if (this->_last_name.empty())
    std::cout << "L'étudiant(e) " << this->_first_name << " a une moyenne de " << this->_average << "." << std::endl;
  else
    std::cout << "L'étudiant(e) " << this->_first_name << " " << this->_last_name << " a une moyenne de " << this->_average << "." << std::endl;
}
// ----------------------------------------------------------

// -------------- Étudiant(e) sérieux(se) -------------------
// Calcul de la moyenne avec la nouvelle note
void SeriousStudent::AddMark(double mark)
{
  this->_average = (this->_average*this->_marks_number + mark)/(this->_marks_number+1.);
  this->_marks_number++;
}
// ----------------------------------------------------------

// -------------- Étudiant(e) en retard ---------------------
// Perd un point à chaque nouvelle note à cause de ses retards !
void LateStudent::AddMark(double mark)
{
  this->_average = (this->_average*this->_marks_number + (mark-1.))/(this->_marks_number+1.);
  this->_marks_number++;
}
// ----------------------------------------------------------

// -------------- Étudiant(e) étourdi(e) --------------------
// Oublie régulièrement de rendre son devoir et n'a donc pas la note méritée !

// Une fonction initialize plus complexe que celle de la classe mère pour fixer this->_dazed
void DazedStudent::Initialize(std::string first_name)
{
  Student::Initialize(first_name);
  // L'étudiant(e) oublie toutes les 4 fois de rendre son devoir
  this->_dazed = 4;
}

// Une fonction initialize plus complexe que celle de la classe mère pour fixer this->_dazed
void DazedStudent::Initialize(std::string first_name, std::string last_name)
{
  Student::Initialize(first_name, last_name);
  // L'étudiant(e) oublie toutes les 4 fois de rendre son devoir
  this->_dazed = 4;
}

// Est-ce que l'étudiant(e) a oublié cette fois ?
int DazedStudent::ForgotOrNot()
{
  if (this->_marks_number % this->_dazed != 0)
    return 1;
  else
    return 0;
}

// Calcul de la moyenne avec sa nouvelle note s'il a pensé à rendre sa copie
void DazedStudent::AddMark(double mark)
{
  if (ForgotOrNot() == 0)
    mark = 0.;

  this->_average = (this->_average*this->_marks_number + mark)/(this->_marks_number+1.);

  this->_marks_number++;
}
// ----------------------------------------------------------
