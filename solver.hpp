#ifndef SOLVER_H
#define SOLVER_H

#include <vector>
#include <assert.h>
#include <iostream>
#include <random>
#include <cmath>

class Solver
{
private:
    std::vector <double> A; //three-diagonal matrix, stored as a linear array: first N + 1 elements form the main diagonal, 
                            //next N elements form the upper subdiagonal, followed by N elements forming the low subdiagonal
    std::vector <double> F; unsigned int N; // N is the number of points
    std::vector <double> p;
    std::vector <double> q;
    double h;
    bool high_order; //if this variable is true, the class uses second-order finite differnce method
                     //if the variable is false, first-order method is used
public: 
    std::vector <double> y;
    Solver (unsigned int _N, std::vector <double> &mesh, double (*P) (double), double (*Q) (double), double (*f) (double),
            double a1, double b1, double g1, double a2, double b2, double g2, bool order) //a1, b1, g1, a2, b2 and g2 are coefficients used in boundary conditions
    {
        N = _N - 1;
        h = mesh[1] - mesh[0];
        assert ("h == 0");
        A = std::vector <double> (3 * N + 1);
        F = std::vector <double> (N + 1);
        if (order)
        {
            A [0] = (Q (mesh [0]) * h * h / 2 - 1) / (h * (1 - P (mesh [0]) * h / 2));
            A [2 * N + 1] = 1 / (h * (1 - P (mesh [0]) * h / 2)); 
            F [0] = g1 + f (mesh [0]) * h / 2 / (1 - h / 2 * P (mesh [0]));
        }
        else
        {
            A [0] = a1 - b1 / h;
            A [2 * N + 1] = b1 / h; 
            F [0] = g1;
        }
        for (unsigned int i = 1; i < N; i++)
        {
            A [N + i] = 1 / h / h - P (mesh [i]) / 2 / h;
            A [i] = Q (mesh [i]) - 2 / h / h;
            A [2 * N + 1 + i] = 1 / h / h + P (mesh [i]) / 2 / h;
            F [i] = f (mesh [i]);
        }
        A [2 * N] = -b2 / h;
        A [N] = a2 + b2 / h;
        F [N] = g2;
        y = std::vector <double> (N + 1);
        p = std::vector <double> (N);
        q = std::vector <double> (N); 
    }
    void Sweep ();
};

void Solver::Sweep ()
{
    p [0] = A [2 * N + 1] / (-A [0]);
    q [0] = F [0] / A [0];
    for (unsigned int i = 0; i < N - 1; i++)
    {
        p [i + 1] = A [2 * N + 1 + i + 1] / (-A [i + 1] - A [N + 1 + i] * p [i]); 
        q [i + 1] = (A [N + 1 + i] * q [i] - F [i + 1]) / (-A [i + 1] - A [N + 1 + i] * p [i]);
    }
    y [N] = (F [N] - A [2 * N] * q [N - 1]) / (A [2 * N] * p [N - 1] + A [N]);
    for (int i = N - 1; i >= 0; i--)
        y [i] = p [i] * y [i + 1] + q [i];
}

class Shooter 
{
private:
    std::vector <double> x;
    std::vector <double> k;
    std::vector <double> q;
    double (*f) (double, double, double);
    double (*g) (double, double, double);
    double y0; //these are initial values 
    double y1; //0 for left end
    double z0; //1 for right end
    double z1;
    double h;
    double eps;
    unsigned int N;
    void Runge_Kutta ();
public:
    std::vector <double> y;
    std::vector <double> z;
    Shooter (unsigned int _N, std::vector <double> &mesh, double (*F) (double, double, double), double (*G) (double, double, double),
                  double _y1, double _z0, double err)
    {
        N = _N;
        h = mesh[1] - mesh[0];
        y1 = _y1;
        z0 = _z0;
        x = mesh;
        y = std::vector <double> (N);
        z = std::vector <double> (N);
        k = std::vector <double> (4);
        q = std::vector <double> (4);
        f = F;
        g = G;
        eps = err;
    }
    void Shoot ();
};

void Shooter::Runge_Kutta ()
{
    y [0] = y0;
    z [0] = z0;
    for (unsigned int i = 0; i < N - 1; i++)
    {
        k [0] = f (x [i], y [i], z [i]);
        q [0] = g (x [i], y [i], z [i]);
        for (unsigned j = 1; j < 3; j++)
        {
            k [j] = f (x [i] + h / 2, y [i] + h * k [j - 1] / 2, z [i] + h * q [j - 1] / 2);
            q [j] = g (x [i] + h / 2, y [i] + h * k [j - 1] / 2, z [i] + h * q [j - 1] / 2);
        }
        k [3] = f (x [i] + h, y [i] + h * k [2], z [i] + h * q [2]);
        q [3] = g (x [i] + h, y [i] + h * k [2], z [i] + h * q [2]);

        y [i + 1] = y [i] + h * (k [0] + 2 * k [1] + 2 * k [2] + k [3]) / 6; 
        z [i + 1] = z [i] + h * (q [0] + 2 * q [1] + 2 * q [2] + q [3]) / 6;
    }
}

void Shooter::Shoot ()
{
    std::random_device rd;
    std::mt19937 gen (rd ());
    double l_bound = -1000.0 + y1;
    double r_bound =  1000.0 + y1;
    std::uniform_real_distribution<> dis (l_bound, r_bound);    
    double tmp1 = dis (gen);
    double tmp2 = 0;
    y0 = tmp1;
    Runge_Kutta ();
    if (y [N - 1] > y1)
        do
        {
            tmp1 = dis (gen);
            y0 = tmp1;
            Runge_Kutta ();
        }
        while (y [N - 1] > y1);
    if (y [N - 1] < y1)
    {
        tmp2 = tmp1;
        do
        {
            tmp2 = dis (gen);
            y0 = tmp2;
            Runge_Kutta ();
        }
        while (y [N - 1] < y1);
    }
    do
    {
        y0 = (tmp1 + tmp2) / 2;
        Runge_Kutta ();            
        if (y [N - 1] < y1)
            tmp1 = y0;
        else 
            tmp2 = y0;
    }
    while (std::abs (y [N - 1] - y1) > eps);
}
#endif
