#include <iostream>
#include <fstream>
#include <chrono>
#include "TimeScheme.h"


int main(int argc, char** argv)
{
   if (argc < 2)
   {
      std::cout << "Please, enter the name of your data file." << std::endl;
      exit(0);
   }
   // ----------------------- Fichier de données --------------------------------
   DataFile* df = new DataFile(argv[1]);

   // ---------------------------- Résolution  ----------------------------------
   // Pointeur contenant toutes les fonctions utiles
   Function* function = new Function(df);
   // Pointeur vers la classe Laplacian (discrétisation en espace)
   Laplacian* lap = new Laplacian(function, df);
   // Pointeur vers la classe TimeScheme (discrétisation en temps=)
   TimeScheme* time_scheme = NULL;

   int int_euler_scheme(1);
   if (df->Get_scheme() == "ImplicitEuler") {int_euler_scheme = 2;}

   switch(int_euler_scheme)
   {
      case 1:
      time_scheme = new EulerScheme(df, lap);
      break;
      case 2:
      time_scheme = new ImplicitEulerScheme(df, lap);
      break;
      default:
      std::cout << "Ce choix n'est pas possible ! Veuillez recommencer !" << std::endl;
      exit(0);
   }

   // ----------------------- Libère la mémoire ---------------------------------
   // Ne pas oublier de détruire tous les pointeurs que nous avons utilisé
   delete df; delete function; delete lap; delete time_scheme;

   return 0;
}
