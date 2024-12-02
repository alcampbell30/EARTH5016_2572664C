%***** RUN 2D MODEL FROM IMAGE ***********************************

% clear workspace
clear all; close all; %clc;

% load model setup from image, interpolate to target grid size
W       = 16e3;     % domain width (must correspond to width of image) [m]
Nx      = 200;      % target grid size z-direction
h       = W/Nx;     % grid spacing based on image width and target grid size
n_units = 9;        % number of rock units contained in image
[units,D,Nz] = ModelFromImage('section.tiff',n_units,W,Nx);

% material properties for each rock unit (update based on your calibration)



matprop = [
% unit conductivity(kT) density(rho)  heat capacity(Cp)  heat production(Hr)
   1	      3.678       2697.6             1000	          4.172         %HE1
   2	      1	          2000	             1000	          1             %Gneiss
   3	      1	          2000	             1000	          1             %Sand
   4	      3.218       2703.5             1000	          5.575         %HE2
   5	      1	          2000	             1000	          1             %Gravel
   6	      1	          2000	             1000	          1             %Clay
   7	      1	          2000	             1000	          1             %Silt
   8	      1	          2000	             1000	          1             %Mud
   9	      1e-6        1000	             1000	          0];           % air/water
         
% get coefficient fields based on spatial distribution of rock units from image
% pay attention if any unit conversion is required!
rho    = reshape(matprop(units,3),Nz,Nx);
Cp     = reshape(matprop(units,4),Nz,Nx);
kT     = reshape(matprop(units,2),Nz,Nx);
Hr     = reshape(matprop(units,5),Nz,Nx);

% continue setting remaining model parameters, then call model routine

%*****  RUN MODEL
run('./ModelFromImage.m');