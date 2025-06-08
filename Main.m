clc
clear
close all

% Wind Turbine Aeroelasticity
% Delft University of Technology
addpath("OurBEM\");

%% 1. Structural Steady Analysis
% Must append NREL5MW.mat size to number of sections in blade section.dat
givenBlade = readtable("Blade\Blade section\Blade section.dat");
R = 63;
r_sections = givenBlade.Radius;

load("NREL5MW.mat","Blade")
Blade.Mass = interp1(Blade.Radius,Blade.Mass,r_sections);
Blade.EIflap = interp1(Blade.Radius,Blade.EIflap,r_sections);
Blade.EIedge = interp1(Blade.Radius,Blade.EIedge,r_sections);
Blade.Twist = givenBlade.AeroTwst;
Blade.Chord = givenBlade.Chord;
Blade.NFoil = givenBlade.AeroNum;
Blade.Radius = r_sections;

dr = givenBlade.DR; 

% Mode shapes
phi_1f = @(r) 0.0622*(r./R).^2 + 1.7254*(r./R).^3 - 3.2452*(r./R).^4 + 4.7131*(r./R).^5 - 2.2555*(r./R).^6;
phi2_1f = @(r) (1/R^2)*(2*0.0622 + 6*1.7254*(r./R) - 12*3.2452*(r./R).^2 + 20*4.7131*(r./R).^3 - 30*2.2555*(r./R).^4);
phi_1e = @(r) 0.3627*(r./R).^2 + 2.5337*(r./R).^3 - 3.5772*(r./R).^4 + 2.376*(r./R).^5 - 0.6952*(r./R).^6;
phi2_1e = @(r) (1/R^2)*(2*0.3627 + 6*2.5337*(r./R) - 12*3.5772*(r./R).^2 + 20*2.376*(r./R).^3 - 30*0.6952*(r./R).^4);

% Mass and stiffness integrals
Mass = Blade.Mass;
Stiffness_Flap = Blade.EIflap;
Stiffness_Edge = Blade.EIedge;
r_struct = Blade.Radius;
dr_struct = dr;

M1f = sum(dr_struct .* Mass .* (phi_1f(r_struct)).^2);
K1f = sum(dr_struct .* Stiffness_Flap .* (phi2_1f(r_struct)).^2);
M1e = sum(dr_struct .* Mass .* (phi_1e(r_struct)).^2);
K1e = sum(dr_struct .* Stiffness_Edge .* (phi2_1e(r_struct)).^2);
damp_ratio = 0.00477465;

M = diag([M1f, M1e])*0.7;
K = diag([K1f, K1e])*0.9;
C = diag([2*damp_ratio*sqrt(M1f*K1f), 2*damp_ratio*sqrt(M1e*K1e)]);

% Compute eigenfrequencies and display
fprintf("Flap freq: %.4f Hz\n", sqrt(K1f/M1f)/(2*pi));
fprintf("Edge freq: %.4f Hz\n", sqrt(K1e/M1e)/(2*pi));


%% 2. Dynamic Inflow
% Initialize inflow conditions
load('STATE');  % Loads WindSpeeds, RtSpeeds, PitchAngles
%import Blade section file
BS = table2array(readtable('Blade/Blade section/Blade section.dat'));
%import Aero data files
Readfiles = dir(fullfile('Blade/Aero data/','*.dat'));
for i=1:length(Readfiles)
    AD{i}=importdata(strcat('Blade/Aero data/',Readfiles(i).name));
end
desired_V = 8;  % for example, 8 m/s

% Find closest matching wind speed
[~, ind] = min(abs(WindSpeeds - desired_V));

% Extract values using that index
Vinf   = WindSpeeds(ind);
Omega  = RtSpeeds(ind) * 2 * pi / 60;  % Convert RPM to rad/s
pitch  = PitchAngles(ind);
twist = deg2rad(Blade.Twist);

dt = 0.005;
t = 0:dt:20;
x = [0;0];
x_dot = [0;0];
x_ddot = [0;0];

for i = 1:length(t)
    
    % Call inplane and out-of-plane velocities
    vout=x_dot(1).*phi_1f(r_struct);
    vin=x_dot(2).*phi_1e(r_struct);

    % Solve BEM
    [Rx, FN, FT, Vind_axial, Vind_tangential] = BEMcode(Vinf,Omega,pitch,vin,vout,BS,AD);

    % Rotate from edge/flap-wise to in/out plane
    ff =  cos(twist) .* FN + sin(twist).* FT;
    fe =  -sin(twist).* FN + cos(twist).* FT;

    % Compute general force and construct F vector
    F1_f=sum(dr_struct .* ff .* phi_1f(r_struct));
    F1_e=sum(dr_struct .* fe .* phi_1e(r_struct));
    F = [F1_f;
         F1_e];
    
    % Compute acceleration using system of equations
    x_ddot = inv(M)*F - inv(M)*C*x_dot - inv(M)*K*x;
    x_dot = x_dot + dt*x_ddot;
    x = x + dt*x_dot;
    xSave(i,:) = x;
    xdotSave(i,:) = x_dot;
end

disp(xSave)
figure;
hold on;
plot(t,xSave(:,1),'-b','LineWidth',1.2,"DisplayName","Out-of-plane displacement");
plot(t,xSave(:,2),'-r','LineWidth',1.2,"DisplayName","In-plane displacement");
grid on;
xlabel("Time [s]");
ylabel("Displacement [m]")
title("Flap and Edge-wise displacement")
legend

figure;
hold on;
plot(t,xdotSave(:,1),'-b','LineWidth',1.2,"DisplayName","Out-of-plane displacement");
plot(t,xdotSave(:,2),'-r','LineWidth',1.2,"DisplayName","In-plane displacement");
grid on;
xlabel("Time [s]");
ylabel("Velocities [m/s]")
title("Flap and Edge-wise velocity")
legend
