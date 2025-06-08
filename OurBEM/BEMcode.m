function [Rx, FN, FT, Vind_axial, Vind_tangential] = BEMcode(Vinf,Omega,Pitch,Vin,Vout,BS,AD)

%Rotor/ Wake Aerodynamics
%Group 28 
%25/03/2025
%Delft University of Technology

%% Airfoil data

%------------------------------------------------
% fixed parameters
%------------------------------------------------
B=3;           %number of blades
R=63;          %rotor radius
hubrad=1.5;    %hub radius
rho=1.225;     %density of air
EPS=0.00001;    %iterative precision tolerance
rRoot = hubrad;

%------------------------------------------------
% Initialization & Iteration
%------------------------------------------------
%initialization: initial value of inductions
a=0;a_prime=0;

NBS=length(BS);    %Number of blade sections
% define vectors for blade section locations and loads in two directions
Rx=zeros(NBS,1);FN=zeros(NBS,1);FT=zeros(NBS,1);

Vind_axial =zeros(NBS,1);
Vind_tangential=zeros(NBS,1);

%% Operational Specs

U0 = Vinf; %[m/s]
omega = Omega; %[rad/s]
rho = 1.225;
q = 0.5 * rho * U0^2;
lambda = omega*R/U0;

%% BEM

% Initialize variables
% CT = zeros(length(lambda),length(mu));
% Cq = zeros(length(lambda),length(mu));
% Cp = zeros(length(lambda),length(mu));
% FAxial = zeros(length(lambda),length(mu));
% FTang = zeros(length(lambda),length(mu));
% aEnd = zeros(length(lambda),length(mu)); 
% aprimeEnd = zeros(length(lambda),length(mu));
% phi_temp = zeros(length(lambda),length(mu));
% alpha_temp = zeros(length(lambda),length(mu));
% Gamma = zeros(length(lambda),length(mu));
% Total_CT = zeros(1,length(lambda));
% Total_CP = zeros(1,length(lambda));

a = 0.3 * ones(1,NBS); % Initialize axial induction factor along blade span
aprime = 0.01 * ones(1,NBS); % Initialize tangential induction factor along blade span
    for j = 1:NBS
    i = 1;
    cond = false;
        while cond == 0

            if i >= 2
            a(i,j) = 0.5*a(i,j) + 0.5*a(i-1,j);
            aprime(i,j) = 0.5*aprime(i,j) + 0.5*aprime(i-1,j);
            else
    
            end
            ADofBS=BS(j,2); % read airfoil number for each section
            r=BS(j,3);      % read radius
            Rx(j)=r;        % record radius
            dr=BS(j,4);     % read segment length
            Theta=BS(j,5);  % read twist angle
            chord=BS(j,6);   % chord length
            Alpha=AD{ADofBS}(:,1);% coefficients table: AOA
            Cl_t=AD{ADofBS}(:,2); % coefficients table: lift coe
            Cd_t=AD{ADofBS}(:,3); % coefficients table: drag coe
            Sigma=chord*B/(2*pi*r); % solidity

            UR = U0*(1-a(i,j))-Vout(j); % Obtain axial velocity at rotor
            UTang = omega*r*(1+aprime(i,j))-Vin(j); % Obtain rotor-induced velocity
            Uapp = sqrt(UR^2 + UTang^2);
    
            phi = atand(UR/UTang); % Inflow angle
            alpha = phi - Theta - Pitch; % AoA
    
            Cl = interp1(Alpha,Cl_t,alpha); % Interpolation to find polars
            Cd = interp1(Alpha,Cd_t,alpha); % Interpolation to find polars
            Cx = Cl*cosd(phi) + Cd*sind(phi);
            Cy = Cl*sind(phi) - Cd*cosd(phi);
            Ct = ((Uapp^2)*Cx*chord*B)/((U0^2)*2*pi*r);
            
            a(i+1,j) = 0.5*(1-sqrt(1-Ct));
    
    %% Glauert Correction
            Ct1 = 1.816;
            if a(i+1,j) > 0.5 % Glauert Correction for induction factors above 0.5
                Ct = 1.816 - 4*(sqrt(1.816)-1)*(1-a(i+1,j));
            else

            end
    
            Ct2 = 2*sqrt(Ct1) - Ct1;
            if Ct >= Ct2
               a(i+1,j) = 1 + (Ct - Ct1)/(4*sqrt(Ct1)-4);           
            else
    
            end
    
            aprime(i+1,j) = ((Uapp^2)*chord*Cy*B*R)/(8*pi*(r^2)*(U0^2)*(1-a(i+1,j))*lambda);
            aDiff = abs(a(i+1,j) - a(i,j));
            aPrimeDiff = abs(aprime(i+1,j) - aprime(i,j));
            cond = (aDiff <= 1e-5 && aPrimeDiff <= 1e-5);
            i = i + 1;
        end
        %% Data storing after convergence
        [a(i,j),aprime(i,j),~] = Prandtl(B,r,R,lambda,a(i,j),aprime(i,j),rRoot);
        UR = U0*(1-a(i,j)); % Obtain axial velocity at rotor
        UTang = omega*r*(1+aprime(i,j)); % Obtain rotor-induced velocity
        Uapp = sqrt(UR^2 + UTang^2);
    
        phiEnd = atand(UR/UTang); % Inflow angle
        alphaEnd = phiEnd - Theta - Pitch; % AoA
        phi_temp(j) = phiEnd;
        alpha_temp(j) = alphaEnd;
    
        Cl = interp1(Alpha,Cl_t,alpha); % Interpolation to find polars
        Cd = interp1(Alpha,Cd_t,alpha); % Interpolation to find polars
        Cx = Cl*cosd(phiEnd) + Cd*sind(phiEnd);
        Cy = Cl*sind(phiEnd) - Cd*cosd(phiEnd);
        CT(j) = 4*a(i,j)*(1-a(i,j));
        % Cq(j) = 4*aprime(i,j)*(1-a(i,j))*lambda*mu(j);
        Cp(j) = 4*a(i,j)*(1-a(i,j))^2;
    
        FAxial(j) = Cx*0.5*rho*(Uapp^2)*chord; % Axial force per newton/m 
        FTang(j) = Cy*0.5*rho*(Uapp^2)*chord; % Tangential force per newton/m
    
        aEnd = a(i,j); 
        aprimeEnd = aprime(i,j);

        % FN(j)=0.5*rho*((r*omega*(1+aprimeEnd))^2+(Vinf*(1-aEnd))^2)*chord*Cx*dr;
        % FT(j)=0.5*rho*((r*omega*(1+aprimeEnd))^2+(Vinf*(1-aEnd))^2)*chord*Cy*dr;
        FN(j) = FAxial(j);
        FT(j) = FTang(j);
        % bending moment
        Mx(j)=FT(j)*r;
    
        % Induced velocities
        Vind_axial(j)      = aEnd.* Vinf;
        Vind_tangential(j) = aprimeEnd .* omega .* Rx(j);
    end
    % TAnnul = FAxial.*dr; %Thrust per annuli of one blade
    % QAnnul = FTang.*dr.*(mu.*R); %Torque per annuli of one blade
    % 
    % QTotal = sum(QAnnul)*B; %Total Torque on rotor (N)
    % TTotal = sum(TAnnul)*B; %Total Thrust on rotor (N)
    % 
    % Total_CT = TTotal/(q*pi*R^2);
    % Total_CP = QTotal*omega/(q*U0*pi*R^2);
end
