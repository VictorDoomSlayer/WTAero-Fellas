function [v0_eff, omega_eff] = compute_effective_inflow(Rx, r, x_dot, phi_f, phi_e, twist, V0, Omega, Vind_axial, Vind_tangential)
% compute_effective_inflow: computes effective axial and tangential velocities
% Inputs:
%   Rx     - BEM blade sections (radius)
%   r      - structural blade sections (radius)
%   x_dot  - [2x1] modal velocities [q1f_dot; q1e_dot]
%   phi_f, phi_e - mode shapes (function handles)
%   twist  - twist angle at each structural location [rad]
%   V0     - freestream wind speed
%   Omega  - rotor rotational speed
%   Vind_axial, Vind_tangential - induced velocities at BEM nodes

% Outputs:
%   v0_eff     - effective axial velocities (out-of-plane)
%   omega_eff  - effective tangential velocities (in-plane)

    % Interpolate induced velocities from Rx (BEM grid) to r (structural grid)
    Vind_axial_interp = interp1(Rx, Vind_axial, r, 'linear', 'extrap');
    Vind_tan_interp   = interp1(Rx, Vind_tangential, r, 'linear', 'extrap');

    % Initialise
    v0_eff = zeros(size(r));
    omega_eff = zeros(size(r));

    % Loop through structural blade sections
    for i = 1:length(r)
        % Modal velocities at position r(i)
        v_flap = x_dot(1) * phi_f(r(i));
        v_edge = x_dot(2) * phi_e(r(i));

        % Rotate from flap/edge to rotor in-plane/out-of-plane
        v_inplane_blade   =  cos(twist(i)) * v_edge - sin(twist(i)) * v_flap;
        v_outplane_blade  =  sin(twist(i)) * v_edge + cos(twist(i)) * v_flap;

        % Compute effective inflow components
        v0_eff(i)     = V0 + Vind_axial_interp(i) - v_outplane_blade;
        omega_eff(i)  = r(i) * Omega + Vind_tan_interp(i) - v_inplane_blade;
    end
end
