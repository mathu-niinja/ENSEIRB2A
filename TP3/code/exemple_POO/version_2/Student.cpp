#include "Student.h"
#include <iostream>

// -------------- Étudiant(e) -------------------------------
// Constructeur par défaut
Student::Student() {}

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

// Ajoute une nouvelle note
void Student::AddMark(double mark)
{
  this->_average = (this->_average*this->_marks_number + mark)/(this->_marks_number+1.);
  this->_marks_number++;
}
// ----------------------------------------------------------
