#include <iostream>

using namespace std;

int main()
{
    int mon_entier(23);
    int *ptr = nullptr; //Un pointeur pouvant contenir l'adresse d'un int
    ptr = &mon_entier; //L'adresse de 'mon_entier' est mise dans le pointeur 'ptr'
    //On dit alors que le pointeur ptr pointe sur mon_entier
    
    cout << "L'adresse de 'mon_entier' est : " << &mon_entier << endl;
    cout << "La valeur de 'ptr' est : " << ptr << endl;
    cout << "La valeur est : " << *ptr << endl;

    int *pointeurInt = nullptr;
    cout << "L'adresse est : " << &pointeurInt << endl;
    cout << "La valeur du pointeur est : " << pointeurInt << endl;


    //On peut faire de même pour n'importe quel type :
    double const *pointeurDoubleConst = nullptr;
    cout << "L'adresse est : " << &pointeurDoubleConst << endl;

    
    return 0;
}