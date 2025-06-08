function Q = compute_generalized_forces(r, dr, FN, FT, phi_f, phi_e, twist)
% Project aerodynamic loads into generalised modal coordinates (Qf, Qe)
    N = length(r);
    Qf = 0;
    Qe = 0;

    for i = 1:N
        % Step 1: Transform from normal/tangential (rotor frame)
        %         to flapwise/edgewise (blade local frame)
        theta = twist(i);  % [rad]

        F_flap =  cos(theta) * FN(i) + sin(theta) * FT(i);
        F_edge = -sin(theta) * FN(i) + cos(theta) * FT(i);

        % Step 2: Project to modal coordinates
        Qf = Qf + F_flap * phi_f(r(i)) * dr(i);
        Qe = Qe + F_edge * phi_e(r(i)) * dr(i);
    end

    Q = [Qf; Qe];
end