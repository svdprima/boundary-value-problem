# Boundary value problem
    
##### MIPT Hometask on Computational Mathematics
        
This programm solves the differential equation y''(x) - y'(x)/x = -2/x^2^ over the domain [0.5; 1] with boundary conditions y'(0.5) = 0, y(1) = 0.
There are three methods implemented and compared:
* order 1 finite difference method
* order 2 finite difference method
* shooting method
                 
The initial value problem in shooting method is solved using the classical 4-order Runge-Kutta method with the following Butcher's table:
                   
0 | 0 | 0 | 0 | 0
--- |--- | --- | --- | --- |
1/2| 1/2 | 0 | 0 | 0
1/2 | 0 | 1/2 | 0 | 0
1 | 0 | 0 | 1 | 0
--- | 1/6 | 1/3 | 1/3 | 1/6

Despite the fact that the programm was written to solve a particular equation, the solver.hpp is still universal for any domain [a;b] and for any two-order differential equation (implying the boundary conditions are y'(a) = y~L~ and y(b) = y~R~) For other boundary conditions the Shooter::Shoot () method should be changed a little.
