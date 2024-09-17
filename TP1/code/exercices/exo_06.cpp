#include <iostream>

using namespace std;

void valeur(int& b)
{
    b = 5;
}

int main()
{
    int a = 3;
    cout << "a avant la fonction valeur " << a << endl;
    valeur(a);
    cout << "a après la fonction valeur " << a << endl;

    return 0;
}