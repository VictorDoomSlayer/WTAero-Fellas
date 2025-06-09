clc
clear
close all

% Wind Turbine Aeroelasticity
% Delft University of Technologyc

addpath("OurBEM/");

% 1. Structural Steady Analysis
% Must append NREL5MW.mat size to number of sections in blade section.dat
givenBlade = readtable("Blade/Blade section/Blade section.dat");
R = 63;
r_sections = givenBlade.Radius;
damp_ratio = 0.00477465;

load("NREL5MW.mat","Blade")
Blade.Mass = interp1(Blade.Radius,Blade.Mass,r_sections);
Blade.EIflap = interp1(Blade.Radius,Blade.EIflap,r_sections);
Blade.EIedge = interp1(Blade.Radius,Blade.EIedge,r_sections);
Blade.Twist = givenBlade.AeroTwst;
Blade.Chord = givenBlade.Chord;
Blade.NFoil = givenBlade.AeroNum;
Blade.Radius = r_sections;

dr=givenBlade.DR;

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


M1f = trapz(r_struct, Mass .* (phi_1f(r_struct)).^2);
K1f = trapz(r_struct, Stiffness_Flap  .* (phi2_1f(r_struct)).^2);
M1e = trapz(r_struct, Mass .* (phi_1e(r_struct)).^2);
K1e = trapz(r_struct, Stiffness_Edge  .* (phi2_1e(r_struct)).^2);

M = diag([M1f, M1e]);
K = diag([K1f, K1e]);
C = diag([2*damp_ratio*sqrt(M1f*K1f), 2*damp_ratio*sqrt(M1e*K1e)]);

% Compute eigenfrequencies and display
fprintf("Flap freq: %.4f Hz\n", sqrt(K1f/M1f)/(2*pi));
fprintf("Edge freq: %.4f Hz\n", sqrt(K1e/M1e)/(2*pi));


%% 2. Dynamic Inflow
% Initialize inflow conditions

 function dxdt = aeroelastic_ode(t, z, M, C, K, Blade, BEMcode, Vinf, Omega, pitch, phi_1f, phi_1e, twist,BS,AD,dr)
    x     = z(1:2);   % Modal displacements [q_f; q_e]
    x_dot = z(3:4);   % Modal velocities [qf_dot; qe_dot]
     
    r_struct = Blade.Radius;

    % Call inplane and out-of-plane velocities
    vout=x_dot(1).*phi_1f(r_struct);
    vin=x_dot(2).*phi_1e(r_struct);

    % Shift to Rotor Axis
    v_axial =  cos(twist + pitch) .* vout - sin(twist + pitch) .* vin;
    v_tang  =  sin(twist + pitch) .* vout + cos(twist + pitch) .* vin;

    % Solve BEM
    [Rx, FN, FT, Vind_axial, Vind_tangential] = BEMcode(Vinf,Omega,rad2deg(pitch),v_tang,v_axial,BS,AD);

    % Rotate from edge/flap-wise to in/out plane
    ff =  cos(twist+pitch) .* FN + sin(twist+pitch).* FT;
    fe =  -sin(twist+pitch).* FN + cos(twist+pitch).* FT;

    % Compute general force and construct F vector
    F1_f=trapz(Blade.Radius, ff./dr .* phi_1f(r_struct));
    F1_e=trapz(Blade.Radius, fe./dr .* phi_1e(r_struct));
    F = [F1_f;
         F1_e];
    
    % Compute acceleration using system of equations
    x_ddot = M\(F - C*x_dot - K*x);
    dxdt = [x_dot; x_ddot];
end

% 3. Solve

%%%%%%%%%%%%%%%%%%%%%%%%%%% Load operational state data %%%%%%%%%%%%%%%%%%%%%%%%%%
load('STATE');  % Loads WindSpeeds, RtSpeeds, PitchAngles
BS = table2array(readtable('Blade/Blade section/Blade section.dat'));
%import Aero data files
Readfiles = dir(fullfile('Blade/Aero data/','*.dat'));
for i=1:length(Readfiles)
    AD{i}=importdata(strcat('Blade/Aero data/',Readfiles(i).name));
end

%%
% Desired wind speed (you can also set this manually)
desired_V = 10;  % for example, 8 m/s

% Find closest matching wind speed
[~, ind] = min(abs(WindSpeeds - desired_V));

% Extract values using that index
Vinf   = WindSpeeds(ind);
Omega  = RtSpeeds(ind) * 2 * pi / 60;  % Convert RPM to rad/s
pitch  = deg2rad(PitchAngles(ind)); % in degrees
twist = deg2rad(Blade.Twist);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z0 = zeros(4, 1);  % No initial displacement or velocity

tspan = [0, 20];  % 0.1 seconds max
odefun = @(t, z) aeroelastic_ode(t, z, M, C, K, Blade, @BEMcode, Vinf, Omega, pitch, phi_1f, phi_1e, twist,BS,AD,dr);
[t, z] = ode45(odefun, tspan, z0);

%
% 4. Postprocess
x1f = z(:,1);      % Flapwise modal displacement
x1e = z(:,2);      % Edgewise modal displacement
dx1f = z(:,3);     % Flapwise velocity
dx1e = z(:,4);     % Edgewise velocity


figure
plot(t, x1f); hold on;
plot(t, x1e); legend('Flapwise', 'Edgewise');
xlabel('Time [s]'); ylabel('Modal displacement [m]');

%%
vin=zeros(17,1);
vout=zeros(17,1);
[Rx, FN, FT, Vind_axial, Vind_tangential] = BEMcode(Vinf,Omega,pitch,vin,vout,BS,AD);

figure(1)
    plot(Rx,FN,'r-o');
    hold on;
    plot(Rx,FT,'b-o');
    hold on
    grid on
    xlabel('Radius(m)');
    ylabel('Loads(N/m)');
    legend('Fn','Ft');

    %%

    %%%%%%%%%%%%%%%%%%%%%%%%%%% Load operational state data %%%%%%%%%%%%%%%%%%%%%%%%%%
load('STATE');  % Loads WindSpeeds, RtSpeeds, PitchAngles
BS = table2array(readtable('Blade/Blade section/Blade section.dat'));
%import Aero data files
Readfiles = dir(fullfile('Blade/Aero data/','*.dat'));
for i=1:length(Readfiles)
    AD{i}=importdata(strcat('Blade/Aero data/',Readfiles(i).name));
end

vels=4:1:24;

tip_def=zeros(length(vels),1);
edge_def=zeros(length(vels),1);

for j=1:length(vels)
% Desired wind speed (you can also set this manually)
desired_V = vels(j);  % for example, 8 m/s

% Find closest matching wind speed
[~, ind] = min(abs(WindSpeeds - desired_V));


% Extract values using that index
Vinf   = WindSpeeds(ind);
Omega  = RtSpeeds(ind) * 2 * pi / 60;  % Convert RPM to rad/s
pitch  = deg2rad(PitchAngles(ind)); % in degrees
twist = deg2rad(Blade.Twist);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z0 = zeros(4, 1);  % No initial displacement or velocity

tspan = [0, 20];  % 0.1 seconds max
odefun = @(t, z) aeroelastic_ode(t, z, M, C, K, Blade, @BEMcode, Vinf, Omega, pitch, phi_1f, phi_1e, twist,BS,AD,dr);
[t, z] = ode45(odefun, tspan, z0);

%
% 4. Postprocess
x1f = z(:,1);      % Flapwise modal displacement
x1e = z(:,2);      % Edgewise modal displacement
dx1f = z(:,3);     % Flapwise velocity
dx1e = z(:,4);     % Edgewise velocity

tip_def(j)=x1f(end);
edge_def(j)=x1e(end);

disp(j)

end

%%
figure
plot(vels, tip_def,'o-'); legend('Flapwise');
grid on
xlabel('Wind Speed [m/s]'); ylabel('Flapwise Tip Deflection [m]');
