clc
clear
close all


% Wind Turbine Aeroelasticity
% Delft University of Technology

% Code for analyzing flapwise and edgewise deformations and dynamics of the
% NREL5MW wind turbine reference

%% 1 Setting the structural parameters

load("NREL5MW.mat")
R = Blade.Radius(end); % Length of Blade
damp_ratio = 0.00477465;


% Following mode factors taken from Slide 32 "Structural Dynamics"
phi_1f = @(r) 0.0622*(r./R).^2 + 1.7254*(r./R).^3 - 3.2452*(r./R).^4 +...
    4.7131*(r./R).^5 - 2.2555*(r./R).^6;
phi2_1f = @(r) (1/R^2)*(6*1.7254*(r./R) - 12*3.2452*(r./R).^2 +...
    20*4.7131*(r./R).^3 - 30*2.2555*(r./R).^4);

phi_1e = @(r) 0.3627*(r./R).^2 + 2.5337*(r./R).^3 - 3.5772*(r./R).^4 +...
    2.376*(r./R).^5 - 0.6952*(r./R).^6;
phi2_1e = @(r) (1/R^2)*(6*2.5337*(r./R) - 12*3.5772*(r./R).^2 +...
    20*2.376*(r./R).^3 - 30*0.6952*(r./R).^4);

% Calculating the eigenfrequency of the FLAP motion

Mass = Blade.Mass;
Stiffness_Flap = Blade.EIflap;
Stiffness_Edge = Blade.EIedge;
M1f = zeros(length(Blade.Radius),1);
K1f = zeros(length(Blade.Radius),1);
M1e = zeros(length(Blade.Radius),1);
K1e = zeros(length(Blade.Radius),1);

dr = diff(Blade.Radius);
dr = [dr; dr(end)];

for i = 1:length(Blade.Radius)
M1f(i) = dr(i)*Mass(i)*(phi_1f(Blade.Radius(i)))^2;
K1f(i) = dr(i)*Stiffness_Flap(i)*(phi2_1f(Blade.Radius(i))^2);
M1e(i) = dr(i)*Mass(i)*(phi_1e(Blade.Radius(i)))^2;
K1e(i) = dr(i)*Stiffness_Edge(i)*(phi2_1e(Blade.Radius(i))^2);
end

M1f = sum(M1f);
K1f = sum(K1f);
M1e = sum(M1e);
K1e = sum(K1e);

K = [K1f 0;
     0  K1e];
M = [M1f 0;
     0  M1e];
C = [2*damp_ratio*sqrt(M1f*K1f)             0;
            0                       2*damp_ratio*sqrt(M1e*K1e)];


omega_flap = sqrt(K1f/M1f);
omega_edge = sqrt(K1e/M1e);

disp("The natural flap frequency is: " + num2str(omega_flap) + " rad/s or " + ...
    num2str(omega_flap/2/pi) + " Hz")
disp("The natural edge frequency is: " + num2str(omega_edge) + " rad/s or " + ...
    num2str(omega_edge/2/pi) + " Hz")

%% 2 Coupling Aerodynamic with Structural Modules

% We need to construct the rest of the code before continuing with question
% 2 of the booklet. We first need to run the BEM with a chosen wind speed
% (TSR) and extract the force per segment (axial and tangential). We then
% must couple it to the structural equation of motion to derive the
% deformation of the section. 

% Each deformation will cause a change of the relative velocity due to a
% flapwise and edgewise displacement. We also need 