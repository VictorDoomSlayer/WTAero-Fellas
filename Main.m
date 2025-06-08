clc
clear
close all

% Wind Turbine Aeroelasticity
% Delft University of Technology

% 1. Structural Parameters
load("NREL5MW.mat")
R = Blade.Radius(end);
damp_ratio = 0.00477465;

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
dr_struct = diff(r_struct); dr_struct(end+1) = dr_struct(end);

M1f = sum(dr_struct .* Mass .* (phi_1f(r_struct)).^2);
K1f = sum(dr_struct .* Stiffness_Flap .* (phi2_1f(r_struct)).^2);
M1e = sum(dr_struct .* Mass .* (phi_1e(r_struct)).^2);
K1e = sum(dr_struct .* Stiffness_Edge .* (phi2_1e(r_struct)).^2);

M = diag([M1f, M1e]);
K = diag([K1f, K1e]);
C = diag([2*damp_ratio*sqrt(M1f*K1f), 2*damp_ratio*sqrt(M1e*K1e)]);

fprintf("Flap freq: %.4f Hz\n", sqrt(K1f/M1f)/(2*pi));
fprintf("Edge freq: %.4f Hz\n", sqrt(K1e/M1e)/(2*pi));


 %%
% 2. Coupling


function dxdt = aeroelastic_ode(t, z, M, C, K, Blade, BEM, Vinf, Omega, pitch, phi_1f, phi_1e, twist)
% ODE function for coupled aeroelastic system
% Converts second-order M*q_ddot + C*q_dot + K*q = F(t) to first-order system


% Split state vector
x     = z(1:2);   % Modal displacements [q_f; q_e]
x_dot = z(3:4);   % Modal velocities [qf_dot; qe_dot]
r_struct = Blade.Radius;

vout=x_dot(1).*phi_1f(r_struct);
vin=x_dot(2).*phi_1e(r_struct);
    
% === STEP 1: Call BEM with steady inflow to get initial induced velocities ===
[Rx, FN, FT, P, Vind_axial, Vind_tangential] = BEM(Vinf, Omega, pitch,vin,vout);
r = Blade.Radius;


% Rotate from flap/edge to rotor in-plane/out-of-plane
ff =  cos(twist) .* FN + sin(twist).* FT;
fe =  -sin(twist).* FN + cos(twist).* FT;

% Rotate from flap/edge to rotor in-plane/out-of-plane
vout_induced =  cos(twist) .* Vind_axial + sin(twist).* Vind_tangential;
vin_induced =  -sin(twist).* Vind_axial + cos(twist).* Vind_tangential;


% Integrate spanwise 
r_struct = Blade.Radius;
dr_struct = diff(r_struct); dr_struct(end+1) = dr_struct(end);
F1_f=sum(dr_struct .* ff .* (phi_1f(r_struct)).^2);
F1_e=sum(dr_struct .* fe .* (phi_1e(r_struct)).^2);


dr=dr_struct;
Q = [F1_f;F1_e]; % returns [Ff; Fe]

% === STEP 5: Solve second-order system and convert to first-order ===
x_ddot = M \ (Q - C * x_dot - K * x);

% === STEP 6: Return dx/dt for ODE solver ===
dxdt = [x_dot; x_ddot];
end



% 3. Solve


%%%%%%%%%%%%%%%%%%%%%%%%%%% Load operational state data %%%%%%%%%%%%%%%%%%%%%%%%%%

load('STATE');  % Loads WindSpeeds, RtSpeeds, PitchAngles
load('NREL5MW.mat');  % Contains Blade.Twist

% Desired wind speed (you can also set this manually)
desired_V = 8;  % for example, 8 m/s

% Find closest matching wind speed
[~, ind] = min(abs(WindSpeeds - desired_V));

% Extract values using that index
Vinf   = WindSpeeds(ind);
Omega  = RtSpeeds(ind) * 2 * pi / 60;  % Convert RPM to rad/s
pitch  = PitchAngles(ind);
twist = deg2rad(Blade.Twist);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%z = [x1f; x1e; dx1f_dt; dx1e_dt];  % 4×1 vector
z0 = zeros(4, 1);  % No initial displacement or velocity


tspan = [0, 0.1];  % 0.1 seconds max
odefun = @(t, z) aeroelastic_ode(t, z, M, C, K, Blade, @BEM, Vinf, Omega, pitch, phi_1f, phi_1e, twist);
%options = odeset('RelTol',1e-3,'AbsTol',1e-4);
[t, z] = ode45(odefun, tspan, z0);


%%
% 4. Postprocess
x1f = z(:,1);      % Flapwise modal displacement
x1e = z(:,2);      % Edgewise modal displacement
dx1f = z(:,3);     % Flapwise velocity
dx1e = z(:,4);     % Edgewise velocity


figure
plot(t, x1f); hold on;
plot(t, x1e); legend('Flapwise', 'Edgewise');
xlabel('Time [s]'); ylabel('Modal displacement [m]');