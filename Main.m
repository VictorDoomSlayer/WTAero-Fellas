clc
clear
close all


% Wind Turbine Aeroelasticity
% Delft University of Technology
%INITIALIZATION Summary of this function goes here

addpath("OurBEM/");

%%%%%%%%%%%%%%%%%%%%%%%%%%% Load operational state data %%%%%%%%%%%%%%%%%%%%%%%%%%
load('STATE');  % Loads WindSpeeds, RtSpeeds, PitchAngles
BS = table2array(readtable('Blade/Blade section/Blade section.dat'));
%import Aero data files
Readfiles = dir(fullfile('Blade/Aero data/','*.dat'));
for i=1:length(Readfiles)
    AD{i}=importdata(strcat('Blade/Aero data/',Readfiles(i).name));
end

% Structural Steady Analysis
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
function [dxdt,F] = aeroelastic_ode(t, z, M, C, K, Blade, BEMcode, Vinf, Omega, pitch, phi_1f, phi_1e, twist,BS,AD,dr,Periodic,Coupling)
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

    if Coupling ==1
        % Shift to Rotor Axis
        v_axial =  cos(twist + pitch) .* vout - sin(twist + pitch) .* vin;
        v_tang  =  sin(twist + pitch) .* vout + cos(twist + pitch) .* vin;
        else
         v_axial =vout*0;
         v_tang=vout*0;
    end


    % Solve BEM
    [Rx, FN, FT, Vind_axial, Vind_tangential] = BEMcode(Vinf,Omega,rad2deg(pitch),v_tang,v_axial,BS,AD);

    % Rotate from edge/flap-wise to in/out plane
    ff =  cos(twist+pitch) .* FN + sin(twist+pitch).* FT;
    fe =  -sin(twist+pitch).* FN + cos(twist+pitch).* FT;

    % Compute general force and construct F vector
    F1_f=trapz(Blade.Radius, ff./dr .* phi_1f(r_struct));
    F1_e=trapz(Blade.Radius, fe./dr .* phi_1e(r_struct));
    F = [F1_f, F1_e];
    
    % Compute acceleration using system of equations
    x_ddot = M\(F' - C*x_dot - K*x);
    dxdt = [x_dot; x_ddot];
 end

%
%% Displacement time-series for V = 15 m/s

Periodic = 0; % 0= Constant inflow, 1= Periodic Wind
Coupling = 1; % aeroelastic coupling activation

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
odefun = @(t, z) aeroelastic_ode(t, z, M, C, K, Blade, @BEMcode, Vinf, Omega, pitch, phi_1f, phi_1e, twist,BS,AD,dr,Periodic,Coupling);
[t, z] = ode45(odefun, tspan, z0);
%[~,F] = cellfun(@(t,z) aeroelastic_ode(t, z.',M, C, K, Blade, @BEMcode, Vinf, Omega, pitch, phi_1f, phi_1e, twist,BS,AD,dr,Periodic),num2cell(t), num2cell(z,2),'uni',0);
%F = cell2mat(F);

% 4. Postprocess
x1f = z(:,1);      % Flapwise modal displacement
x1e = z(:,2);      % Edgewise modal displacement
dx1f = z(:,3);     % Flapwise velocity
dx1e = z(:,4);     % Edgewise velocity

% 5. Compute root bending moments
Moments = zeros(2, length(t));
for i = 1:length(t)
    Moments(:,i) = compute_moments(z(i,1:2)', z(i,3:4)', Blade, phi_1f, phi_1e, pitch, twist, Vinf, Omega, BS, AD, dr);
end

% Sanity check for BEM
vin=zeros(17,1);
vout=zeros(17,1);
[Rx, FN, FT, Vind_axial, Vind_tangential] = BEMcode(Vinf,Omega,pitch,vin,vout,BS,AD);


% Plot deflections (steady)
figure(1)
subplot(2,1,1)
plot(t, x1f); grid on;
title('Tip deflection')
ylabel('Modal displacement [m]');

subplot(2,1,2)
plot(t, x1e); 
title('Edge deflection')
xlabel('Time [s]'); ylabel('Modal displacement [m]');

%% Out-of-plane tip deflection from 4 m/s to 24 m/s

vels=4:1:25;
Periodic = 0;
Coupling =1; 
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
odefun = @(t, z) aeroelastic_ode(t, z, M, C, K, Blade, @BEMcode, Vinf, Omega, pitch, phi_1f, phi_1e, twist,BS,AD,dr,Periodic,Coupling);
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

%% Plot
figure

plot(vels, tip_def,'o-',Color='b',LineWidth=0.8);
grid on

xlabel('Wind Speed [m/s]'); ylabel('Flapwise Tip Deflection [m]');
xlim([3 25.5])
%% Periodic wind inflow

Periodic = 1; % 0= Constant inflow, 1= Periodic Wind
Coupling = 1; % aeroelastic coupling activation


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
clear z 
tspan = [0, 50];  % 0.1 seconds max
odefun = @(t, z) aeroelastic_ode(t, z, M, C, K, Blade, @BEMcode, Vinf, Omega, pitch, phi_1f, phi_1e, twist,BS,AD,dr,Periodic,Coupling);
[t, z] = ode45(odefun, tspan, z0);


% Postprocess
x1f = z(:,1);      % Flapwise modal displacement
x1e = z(:,2);      % Edgewise modal displacement
dx1f = z(:,3);     % Flapwise velocity
dx1e = z(:,4);     % Edgewise velocity

% 5. Compute root bending moments
Moments = zeros(2, length(t));
for i = 1:length(t)
    Moments(:,i) = compute_moments(z(i,1:2)', z(i,3:4)', Blade, phi_1f, phi_1e, pitch, twist, Vinf, Omega, BS, AD, dr);
end

%%
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
%%%%%%%%%%%%%%% plot root bending moments %%%%%%%%%%%%
figure;
subplot(2,1,1)
plot(t, M(1,:),'LineWidth',1.2);
 ylabel('Moment [Nm]');
grid on;
title("Flapwise root bending moment")

subplot(2,1,2)
plot(t, M(2,:),'LineWidth',1.2);
xlabel('Time [s]');   ylabel('Moment [Nm]');
grid on;
title("Edgewise root bending moment")

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

%titles={'Flapwise root bending moment','Edgewise root bending moment'};

signal=[M(1,:) ; M(2,:)];
for i=1:2

% Define uniform time vector
t_uniform = linspace(t(1), t(end), length(t));

% Interpolate signal
y_uniform = interp1(t, signal(i,:), t_uniform, 'pchip');
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
%subplot(2,1,i)
semilogy(f,P1) 
grid on
hold on
%title(titles{i})
xlabel("f (Hz)")
ylabel("|P1(f)|")
xlim([0 3])

%xline(0.7022, '--b', 'f_e (1.079 Hz)', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'right', 'FontSize', 7);

end
xline(1.079, '--b', 'f_e (1.079 Hz)', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'right', 'FontSize', 7);
xline(1.267/2/pi, 'k--', '1P','LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'center','FontSize', 10);
xline(2.534/2/pi, 'k--', '2P','LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'bottom','LabelHorizontalAlignment', 'center', 'FontSize', 10);
xline(3.801/2/pi, 'k--', '3P','LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'center', 'FontSize', 10);
legend('Flapwise','Edgewise')
title('Frequency response, root bending moment')

%% Coupling with Periodic wind inflow
clc

for j=0:1
Periodic = 1; % 0= Constant inflow, 1= Periodic Wind
Coupling = j; % aeroelastic coupling activation


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
odefun = @(t, z) aeroelastic_ode(t, z, M, C, K, Blade, @BEMcode, Vinf, Omega, pitch, phi_1f, phi_1e, twist,BS,AD,dr,Periodic,Coupling);
[t, z] = ode45(odefun, tspan, z0);


% Postprocess
x1f = z(:,1);      % Flapwise modal displacement
x1e = z(:,2);      % Edgewise modal displacement
dx1f = z(:,3);     % Flapwise velocity
dx1e = z(:,4);     % Edgewise velocity

% 5. Compute root bending moments
Moments = zeros(2, length(t));
for i = 1:length(t)
    Moments(:,i) = compute_moments(z(i,1:2)', z(i,3:4)', Blade, phi_1f, phi_1e, pitch, twist, Vinf, Omega, BS, AD, dr);
end


%%%%%%%%%%%%%%% plot displacement %%%%%%%%%%%%
figure(1);
plot(t, x1f,'LineWidth',1.2);
 ylabel('Displacement [m]');
grid on;
title("Flapwise response")
hold on

%%%%%%%%%%%%%%% plot moments %%%%%%%%%%%%
figure(2);
subplot(2,1,1)
plot(t, Moments(1,:),'LineWidth',1.2);
 ylabel('Moment [Nm]');
grid on;
hold on
title("Flapwise root bending moment")

subplot(2,1,2)
plot(t, Moments(2,:),'LineWidth',1.2);
xlabel('Time [s]');   ylabel('Moment [Nm]');
grid on;
hold on
title("Edgewise root bending moment")


%%%%%%%%%%%%%%% plot frequency response%%%%%%%%%%%%
figure(3);
titles={'Flapwise root bending moment','Edgewise root bending moment'};
for k = 1:2 
    signal=[Moments(1,:) ; Moments(2,:)];
    % Define uniform time vector
    t_uniform = linspace(t(1), t(end), length(t));
    % Interpolate signal
    y_uniform = interp1(t, signal(k,:), t_uniform, 'pchip');
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
    subplot(2,1,k)
    semilogy(f,P1) 
    grid on
    hold on
    title(titles{k})
    xlabel("f (Hz)")
    ylabel("|P1(f)|")
    xlim([0 4])
    
end
end

figure(3)
subplot(2,1,2)
xline(1.079, '--b', 'f_e (1.079 Hz)', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'right', 'FontSize', 7);
xline(1.267/2/pi, 'k--', '1P','LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'center','FontSize', 10);
xline(2.534/2/pi, 'k--', '2P','LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'bottom','LabelHorizontalAlignment', 'center', 'FontSize', 10);
xline(3.801/2/pi, 'k--', '3P','LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'center', 'FontSize', 10);
legend('Without Coupling','Coupled')