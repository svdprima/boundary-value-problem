#include <vector>
#include <cmath>
#include <fstream>
#include "solver.hpp"

double solution (double x) //f is the exact solution
{
    return -2 * x * x + 2 + log (x);
}

double p (double x) // p and q are variable coeffs
{
    return -1 / x;
}

double q (double x)
{
    return 0;
}

double f (double x)
{
    return -2 / x / x;
}

double F (double x, double u, double v) // F and G are functions in a system of two first-order diff equations
{
    return v;
}

double G (double x, double u, double v)
{
    return v / x - 2 / x / x;
}

int main ()
{
    const unsigned int N = 10;
    const double a = 0.5;
    const double b = 1.0;
    const double h = (b - a) / (N - 1);
    std::vector <double> mesh = std::vector <double> (N);
    for (unsigned int i = 0; i < N; i++)
        mesh[i] = a + h * i;
    {
        std::ofstream output;
        output.open ("exact.txt");
        for (unsigned int i = 0; i < N; i++)
            output << mesh [i] << " " << solution (mesh [i]) << std::endl;
        output.close ();
    }
    Solver S1 = Solver (N, mesh, p, q, f, 0, 1, 0, 1, 0, 0, 0);
    S1.Sweep ();
    {
        std::ofstream output;
        output.open ("first_order.txt");
        for (unsigned int i = 0; i < N; i++)
            output << mesh [i] << " " << S1.y [i] << std::endl;
        output.close ();
    }
    Solver S2 = Solver (N, mesh, p, q, f, 0, 1, 0, 1, 0, 0, 1);
    S2.Sweep ();
    {
        std::ofstream output;
        output.open ("second_order.txt");
        for (unsigned int i = 0; i < N; i++)
            output << mesh [i] << " " << S2.y [i] << std::endl;
        output.close ();
    }
    double error = 0.00001;
    Shooter S3 = Shooter (N, mesh, F, G, 0.0, 0.0, error);
    S3.Shoot ();
    {
        std::ofstream output;
        output.open ("shooting.txt");
        for (unsigned int i = 0; i < N; i++)
            output << mesh [i] << " " << S3.y [i] << std::endl;
            output.close ();
    }
    return 0;
}
