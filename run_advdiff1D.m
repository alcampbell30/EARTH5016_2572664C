%***** RUN 1D ADVECTION DIFFUSION MODEL ***********************************

% clear workspace
clear all; close all; %clc;

% set model parameters
W  = 500;          % domain width [m]
w  = 50;           % magma width [m]
N  = 200;          % grid size
dx = W/N;          % grid spacing

T0   = 100;        % initial background temperature [C]
dT   = 1000;       % initial temperature peak amplitude [C]
wT   = 20;         % initial temperature peak width [m]
k0   = 0e-6;       % heat diffusivity [m2/s]
u0   = 1e-6;       % advection speed [m/s]
BC   = 'periodic'; % boundary condition option flag ('insulating', 'periodic')
ADVN = 'WENO5';    % advection scheme ('UPW1', 'CFD2', 'UPW3', 'WENO5')
TINT = 'FE1';      % time integration scheme ('FE1','RK2','RK4')

yr    = 3600*24*365;  % seconds per year [s]
tend  = W/u0;         % stopping time [s]
CFL   = 1/16;         % Time step limiter
nop   = 100;          % output figure produced every 'nop' steps

%*****  RUN MODEL
run('./advdiff1D_solution.m');