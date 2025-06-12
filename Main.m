clc
clear
close all
% sebas
% Wind Turbine Aeroelasticity
% Delft University of Technologyc

addpath("OurBEM/");

%% Structural Steady Analysis
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
omega_flap = sqrt(K1f/M1f)/(2*pi); %[hz]
omega_edge = sqrt(K1e/M1e)/(2*pi); %[hz]
fprintf("Flap freq: %.4f Hz\n", omega_flap);
fprintf("Edge freq: %.4f Hz\n", omega_edge);


%% Dynamic Inflow
% Setup aeroelastic equation of motion
 function dxdt = aeroelastic_ode(t, z, M, C, K, Blade, BEMcode, Vinf, Omega, pitch, phi_1f, phi_1e, twist,BS,AD,dr,Periodic)
    x     = z(1:2);   % Modal displacements [q_f; q_e]
    x_dot = z(3:4);   % Modal velocities [qf_dot; qe_dot]
     
    r_struct = Blade.Radius;

    if Periodic == 1
        Vinf = 15 + 0.5*cos(1.267*t) + 0.085*cos(2.534*t) + 0.015*cos(3.801*t);
    else

    end

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

%%%%%%%%%%%%%%%%%%%%%%%%%%% Load operational state data %%%%%%%%%%%%%%%%%%%%%%%%%%
load('STATE');  % Loads WindSpeeds, RtSpeeds, PitchAngles
BS = table2array(readtable('Blade/Blade section/Blade section.dat'));
%import Aero data files
Readfiles = dir(fullfile('Blade/Aero data/','*.dat'));
for i=1:length(Readfiles)
    AD{i}=importdata(strcat('Blade/Aero data/',Readfiles(i).name));
end

%% Displacement time-series for V = 15 m/s
Periodic = 0; % 0= Constant inflow, 1= Periodic Wind
% Desired wind speed (you can also set this manually)
desired_V = 15;  % for example, 8 m/s

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
odefun = @(t, z) aeroelastic_ode(t, z, M, C, K, Blade, @BEMcode, Vinf, Omega, pitch, phi_1f, phi_1e, twist,BS,AD,dr,Periodic);
[t, z] = ode45(odefun, tspan, z0);

% 4. Postprocess
x1f = z(:,1);      % Flapwise modal displacement
x1e = z(:,2);      % Edgewise modal displacement
dx1f = z(:,3);     % Flapwise velocity
dx1e = z(:,4);     % Edgewise velocity


% Sanity check for BEM
vin=zeros(17,1);
vout=zeros(17,1);
[Rx, FN, FT, Vind_axial, Vind_tangential] = BEMcode(Vinf,Omega,pitch,vin,vout,BS,AD);


% Plot deflections (steady)
figure(1)
subplot(2,1,1)
plot(t, x1f); grid on;
ylim([0 3.5])
title('Tip deflection')
ylabel('Modal displacement [m]');

subplot(2,1,2)
plot(t, x1e); 
title('Edge deflection')
xlabel('Time [s]'); ylabel('Modal displacement [m]');

%% Out-of-plane tip deflection from 4 m/s to 24 m/s

vels=4:1:24;
Periodic = 0;
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
odefun = @(t, z) aeroelastic_ode(t, z, M, C, K, Blade, @BEMcode, Vinf, Omega, pitch, phi_1f, phi_1e, twist,BS,AD,dr,Periodic);
[t, z] = ode45(odefun, tspan, z0);


% 4. Postprocess
x1f = z(:,1);      % Flapwise modal displacement
x1e = z(:,2);      % Edgewise modal displacement
dx1f = z(:,3);     % Flapwise velocity
dx1e = z(:,4);     % Edgewise velocity

tip_def(j)=x1f(end);
edge_def(j)=x1e(end);

disp(j)

end

figure
plot(vels, tip_def,'o-'); legend('Flapwise');
grid on
xlabel('Wind Speed [m/s]'); ylabel('Flapwise Tip Deflection [m]');

%% Periodic wind inflow

Periodic = 1; % 0= Constant inflow, 1= Periodic Wind

% Desired wind speed (you can also set this manually)
desired_V = 15;

% Find closest matching wind speed
[~, ind] = min(abs(WindSpeeds - desired_V));

% Extract values using that index
Vinf   = WindSpeeds(ind);
Omega  = RtSpeeds(ind) * 2 * pi / 60;  % Convert RPM to rad/s
pitch  = deg2rad(PitchAngles(ind)); % in degrees
twist = deg2rad(Blade.Twist);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z0 = zeros(4, 1);  % No initial displacement or velocity

tspan = [0, 50];  % 0.1 seconds max
odefun = @(t, z) aeroelastic_ode(t, z, M, C, K, Blade, @BEMcode, Vinf, Omega, pitch, phi_1f, phi_1e, twist,BS,AD,dr,Periodic);
[t, z] = ode45(odefun, tspan, z0);

% Postprocess
x1f = z(:,1);      % Flapwise modal displacement
x1e = z(:,2);      % Edgewise modal displacement
dx1f = z(:,3);     % Flapwise velocity
dx1e = z(:,4);     % Edgewise velocity

%%%%%%%%%%%%%%% plot displacement %%%%%%%%%%%%
figure;
subplot(2,1,1)
plot(t, x1f,'LineWidth',1.2);
 ylabel('Displacement [m]');
grid on;
title("Flapwise response")

subplot(2,1,2)
plot(t, x1e,'LineWidth',1.2);
xlabel('Time [s]');  ylabel('Displacement [m]');
grid on;
title("Edgewise response")

%%
figure;
subplot(2,1,1)
plot(t, dx1f,'LineWidth',1.2); hold on;
ylabel('Modal velocity [m/s]');
grid on
title("Flapwise velocity response")
subplot(2,1,2)
plot(t, dx1e,'LineWidth',1.2);
xlabel('Time [s]'); ylabel('Modal velocity [m/s]');
grid on;
title("Edgewise velocity response")
%%
plotFrequencySpectra(t, x1f, x1e, dx1f, dx1e, Omega, omega_flap, omega_edge);
%% New plotting routine for Freq response

for i=1:2
signal=[x1f , x1e];
% Define uniform time vector
t_uniform = linspace(t(1), t(end), length(t));

% Interpolate signal
y_uniform = interp1(t, signal(:,i), t_uniform, 'pchip');
%y_uniform = y_uniform - mean(y_uniform);

dt=mean(diff(t_uniform));% approximate mean timestep
L = length(t_uniform);       % Number of points
Fs = 1/dt;                   % Sampling frequency
window = hann(L)';
Y = fft(y_uniform .* window);% Compute FFT}


%Y = fft(y_uniform);
P2 = abs(Y/L);                % Magnitude spectrum (normalised)

P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs/L*(0:(L/2));

%figure(2)
subplot(2,1,i)
semilogy(f,P1) 
grid on
hold on
title("Bending root moment Frequency Response")
xlabel("f (Hz)")
ylabel("|P1(f)|")
xlim([0 5])

%xline(0.7022, '--b', 'f_e (1.079 Hz)', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'right', 'FontSize', 7);
xline(1.079, '--k', 'f_e (1.079 Hz)', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'right', 'FontSize', 7);
end

legend('Flapwise', 'Edgewise')