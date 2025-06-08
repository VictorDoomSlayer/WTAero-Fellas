function [Rx, FN, FT, P, Vind_axial, Vind_tangential] = BEM_vector(v0, omega, pitch)
%------------------------------------------------
% BEM_vector - Vectorised BEM tailored for aeroelastic coupling
% TU Delft - Wind Turbine Aeroelasticity
%------------------------------------------------

% Constants
B       = 3;           % number of blades
R       = 63;          % rotor radius [m]
hubrad  = 1.5;         % hub radius [m]
rho     = 1.225;       % air density [kg/m^3]
EPS     = 1e-5;        % convergence tolerance

% Load blade section data
BS = table2array(readtable('Blade/Blade section/Expanded_Blade_Section_Table__48_points_.csv'));
Readfiles = dir(fullfile('Blade/Aero data/', '*.dat'));
for j = 1:length(Readfiles)
    AD{j} = importdata(fullfile('Blade/Aero data/', Readfiles(j).name));
end

% Discretisation
N = size(BS, 1);
Rx = zeros(N,1); FN = zeros(N,1); FT = zeros(N,1); Mx = zeros(N,1);
a_out = zeros(N,1); a_prime_out = zeros(N,1);

% Broadcast scalars if needed
if isscalar(v0),     v0 = v0 * ones(N,1); end
if isscalar(omega),  omega = omega * ones(N,1); end

% Loop over blade elements
for i = 1:N
    r       = BS(i,2);     Rx(i) = r;
    dr      = BS(i,4);
    theta   = BS(i,5);     % twist angle in deg
    chord   = BS(i,6);
    af_idx  = BS(i,3);     % airfoil index

    % Airfoil data
    alpha_tab = AD{af_idx}(:,1);
    cl_tab    = AD{af_idx}(:,2);
    cd_tab    = AD{af_idx}(:,3);

    sigma = chord * B / (2*pi*r);
    a = 0.3; a_prime = 0.01;

    for iter = 1:100
        omega_i = max(omega(i), 1e-3);  % avoid div-by-zero
        phi = atan((1-a)*v0(i) / ((1+a_prime)*r*omega_i));
        phi_deg = rad2deg(phi);
        alpha = real(phi_deg - theta - pitch);
        if isnan(alpha), alpha = 0; end

        cl = interp1(alpha_tab, cl_tab, alpha, 'linear', 'extrap');
        cd = interp1(alpha_tab, cd_tab, alpha, 'linear', 'extrap');

        cn = cl * cosd(phi_deg) + cd * sind(phi_deg);
        ct = cl * sind(phi_deg) - cd * cosd(phi_deg);

        % Prandtl loss
        f_tip = (B/2)*(R - r)/(r * sind(phi_deg));
        f_hub = (B/2)*(r - hubrad)/(r * sind(phi_deg));
        F_tip = (2/pi) * acos(exp(-f_tip));
        F_hub = (2/pi) * acos(exp(-f_hub));
        F = max(F_tip * F_hub, 1e-4);

        % Glauert correction
        ac = 0.2;
        if a > ac
            K = 4*F*sind(phi_deg)^2 / (sigma*cn);
            a_new = 0.5*(2 + K*(1 - 2*ac) - sqrt((K*(1 - 2*ac) + 2)^2 + 4*(K*ac^2 - 1)));
        else
            a_new = 1 / (4*F*sind(phi_deg)^2 / (sigma*cn) + 1);
        end

        a_prime_new = 1 / (4*F*sind(phi_deg)*cosd(phi_deg) / (sigma*ct) - 1);

        if abs(a_new - a) < EPS && abs(a_prime_new - a_prime) < EPS
            break;
        end
        a = a_new; a_prime = a_prime_new;
    end

    a_out(i) = a;
    a_prime_out(i) = a_prime;

    Vrel2 = (v0(i)*(1-a))^2 + (r*omega(i)*(1+a_prime))^2;
    FN(i) = 0.5 * rho * Vrel2 * chord * cn * dr;
    FT(i) = 0.5 * rho * Vrel2 * chord * ct * dr;
    Mx(i) = FT(i) * r;
end

% Total power output
P = sum(Mx .* omega) * 3 * 0.944;

% Induced velocities
Vind_axial      = a_out       .* v0;
Vind_tangential = a_prime_out .* omega .* Rx;

end
