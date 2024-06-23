clc;
clear all;
close all;

U = 5;
alpha = [1; 2; 5; 3; 4];
j = 2;
N_0 = 10.^(-20.4);
del = 1;
d_0 = 24.5e6;
% sigma = 1;
t_0 = 0.01

Bmin = 1e6;
Bmax = 30e6;
Pmin = 0.001;
Pmax = 1;
tmin = 0.001;
tmax = 0.010;

[B,P] = fn(U,Bmin,Bmax,Pmin,Pmax,t_0,tmin,tmax,alpha,j,N_0,del,d_0);
