#include "Student.h"
#include <string>
#include <iostream>
#include <vector>

int main()
{
	// Un groupe de 4 étudiants
	std::vector<Student*> students(4);

	students[0] = new Student(); // Premier(e) étudiant(e)
	students[0]->Initialize("Gabriel", "Guerin"); // Son nom de famille et son prénom sont connus

	students[1] = new Student(); // Second(e) étudiant(e)
	students[1]->Initialize("Jules", ""); // Son prénom slmt est connu

	students[2] = new Student(); // Troisième étudiant(e)
	students[2]->Initialize("Louise", ""); // Son prénom slmt est connu

	students[3] = new Student(); // Quatrième étudiant(e)
	students[3]->Initialize("Alice", "Perrin"); // Son nom de famille et son prénom sont connus

	// Toutes les notes qu'ils ont obtenues
	std::vector<double> marks = {12.3, 14.5, 16.8, 9.5, 10.3, 18.1, 14.3, 16.1, 19.8, 15};

	for (int i = 0 ; i < students.size() ; ++i)
	{
		for (int j = 0 ; j < marks.size() ; ++j)
		{
			students[i]->AddMark(marks[j]);
		}
		students[i]->Print();
	}


	for (int i = 0 ; i < students.size() ; ++i)
	{
		delete students[i];
	}

	return 0;
}