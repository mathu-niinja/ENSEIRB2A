#include <cmath>
#include <iostream>

using namespace std;

void euler(double &y, double tn, double dt, double (&f)(double, double))
{
    y += dt * f(tn, y);
}

double f_ty2(double t, double y)
{
    return t * y * y;
}

int main()
{
    cout.precision(15);
    double dt(0.001), Tfinal(1.), y0(-2.), y(y0), tn;
    int N = int(round(Tfinal / dt));

    for (int i = 1; i < N + 1; i++)
    {
        tn = i * dt;
        euler(y, tn, dt, f_ty2);
    }

    double y_ext_T_final = 2. * y0 / (2. - y0 * Tfinal * Tfinal);
    
    cout << "y'=t*y^2 avec y(0)=" << y0 << endl;
    cout << "y(Tf)=" << y << " and y_exacte(Tf)=" << y_ext_T_final << endl;

    return 0;
}