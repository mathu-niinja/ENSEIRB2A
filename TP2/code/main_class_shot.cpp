#include <fstream>
#include "Shot.h"

using namespace std;

int main()
{
    // Initial conditions
    double x0 = 0;
    double y0 = 1.8;
    double normv0 = 7.2;
    double anglev0 = 0.6;

    // Build shot
    double mass = 4000; // g
    Shot shot(mass, y0, normv0, anglev0);
    
    // Throw shot
    shot.computeTrajectory();
    // Print shot after
    cout << shot << endl;
    // Write in a file
    string name_file("Results/solution_1_shot.txt");
    shot.saveSolution(name_file);

    return 0;
}