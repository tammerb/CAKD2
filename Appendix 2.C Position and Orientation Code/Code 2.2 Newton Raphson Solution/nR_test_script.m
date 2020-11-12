% Clear Everything
clearvars;          % Clear all variables from the workspace
clc;                % Clear Command Window

%This case use a 3 dimension as an example (This code works for n-dimensions)

%define variable
syms x1 x2 x3;      % Define varible base on number of dimension (x1 ... xn)
sym_x = [x1 x2 x3].';% System of variable x

%1. Test a system of linear equation
disp('Test 1')
xtol = 0.001;  % Enter error tolerance in satisfying F(x)=0
x0 = [0.1 0.1 0.1].';   % Enter initial solution estimate 1xn matrix
sym_F = [8*x1 4*x2 x3-4].';   % Enter system of function of F(x)
x = newtonRaphson(sym_x, xtol, x0, sym_F)

%2. Test a system of nonlinear equation
disp('Test 2')
xtol = 0.001;  % Enter error tolerance in satisfying F(x)=0
x0 = [0.1 0.1 0.1].';   % Enter initial solution estimate 1xn matrix
sym_F = [cos(x1)*sin(x2) cos(x2)^2 4*x1*cos(x3)].';   % Enter system of function of F(x)
x = newtonRaphson(sym_x, xtol, x0, sym_F)

%3. Test a system of 3rd oreder polynomials
disp('Test 3')
xtol = 0.001;  % Enter error tolerance in satisfying F(x)=0
x0 = [0.1 0.1 0.1].';   % Enter initial solution estimate 1xn matrix
sym_F = [2*x1^2+x2-x3 x3^3-4 x2^2+4*x2-1].';   % Enter system of function of F(x)
x = newtonRaphson(sym_x, xtol, x0, sym_F)

%4 Test the algorithm with a trig function that doesn't converge
disp('Test 4')
xtol = 0.001;  % Enter error tolerance in satisfying F(x)=0
x0 = [0.1 0.1 0.1].';   % Enter initial solution estimate 1xn matrix
sym_F = [cos(x1)*tan(x2) (cos(x2))^2 4*x1*cos(x3)].';   % Enter system of function of F(x)
x = newtonRaphson(sym_x, xtol, x0, sym_F)