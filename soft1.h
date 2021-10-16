#include<stdio.h>


int i, j, k, count;
int m = 500, n = 100, n1, n2, m1, m2, m3, m4, m5, m6;


double x[1001][501], y[1001][501], fakeomega[1001][501], omega[1001][501], omegaold[1001][501], u[1001][501], v[1001][501];
double stream[1001][501], sy[1001][501], syold[1001][501], fakesy[1001][501];
double phi[1001][501], phiold[1001][501], fakephi[1001][501];
double f[1001][501], f1[1001][501], f2[1001][501];
double sp[1001][501], spold[1001][501], fakesp[1001][501];

double flag, app_error, error, sum, dx, dy, flag, dif;

double wphi = 1.5;
double wsp = 0.1;
double wsy = 0.9;
double womega = 0.8;


double uin = 0.25;
double re = 5.0e-3;
double sc = 1.0e5;
double kappa = 16.0;
double d = 0.2;
double alpha = 2.0;


double app_pot;
double kappap = 10.0; 
double length = 5.0;
double height = 1.0;
