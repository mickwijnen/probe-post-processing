function [iIth, DL, xi] = paraBRL(V, iVar, eVar, const, probe)
% Generates theoretical BRL ion-current according to Chen's parametrization
%
% Based on paper; Title:    Langmuir probe analysis for high density plasma
%                 Author:   Francis F. Chen
%                 Journal:  Physics of Plasmas                  
%                 Volume:   8
%                 Number:   6
%                 Year:     2001
%                 DOI:      10.1063/1.1368874
%
% All equations and table numbers refer to the above paper.
%
% Output:           iIth: theoretical ion current according to BRL theory
%                   DL:   Debye length [mm]
%                   xi:   ratio of probe radius and Debye length
% Input:
%                   V:      voltage data array
%                   iVar:   ion variables, n iVs
%                   eVar:   electron variables, Te eVs
%                   values: struct with values for (Te, n, iVs, eVs)
%                   const:  struct with (natural) constants, e
%                   probe:  struct with probe dimensions, Rp Lp S AR
%                   


%------------------------Constants & Variables---------------------------%

Te  = eVar(1);              % [eV] electron temperature
eVs = eVar(2);              % [V] electron plasma potential

n   = iVar(1)*1e17;         % [m^-3] plasma density
iVs = iVar(2);              % [V] ion plasma potential

e    = const.e;             % [C] elementary charge
eps0 = const.eps0;          % [F/m] vacuum permittivity
mi   = const.mi;   

Rp = probe.Rp;              % [mm] probe radius
S  = probe.S;               % [mm^2] probe surface

Va = sqrt(e*Te/2/pi/mi);    % thermal ion velocity

%--------------------Calculate auxilliary variables----------------------%

DL = 1000*sqrt(eps0*Te/e/n); % [mm] Debye length in 

xi = Rp/DL;                  % probe radius to Debye length ratio

Jr = e*S*Va;                 % [Cm^3/s] random ion current per density
    
eta = (iVs-V)/Te;            % non-dimensional potential

eta = eta.*(eta > 0);        % make all negative values 0

%--------------------Generate parametric coefficients--------------------%

if xi > 3
    
    % equation(9)
    funA = @(q) q(1)+q(2)*(xi-q(3))^q(4)*exp(-q(5)*(xi-q(3))^q(6));
    funB = funA;
    funC = @(q) q(1)+q(2)*exp(-q(3)*log(xi-q(4)))+q(5)*(1-q(6)*log(xi));
    funD = funA;
    
    % Table II
    coef = [ 1.142 19.027 3.000 1.433 4.164 0.252;
             0.530  0.970 3.000 1.110 2.120 0.350;
             0.000  1.000 3.000 1.950 1.270 0.035;
             0.000  2.650 2.960 0.376 1.940 0.234; ];
    
    A = funA(coef(1,:));
    B = funB(coef(2,:));
    C = funC(coef(3,:));
    D = funD(coef(4,:));
       
else
    
    % equation(10)
    funA=@(q) q(1)+1/(1/q(2)/xi^q(3)-1/q(4)/log(xi/q(5)));
    funB=@(q) q(1)+q(2)*xi^q(3)*exp(-q(4)*xi^q(5));
    funC=@(q) q(1)+q(2)*xi^-q(3);
    funD=funB;
  
    % Table IV
    coef = [ 1.12 0.00034 6.87 0.145 110.0;
             0.50 0.00080 1.50 0.180 0.800;
             1.07 0.95000 1.01   NaN   NaN;
             0.05 1.54000 0.30 1.135 0.370; ];
    
    A = funA(coef(1,:));
    B = funB(coef(2,:));
    C = funC(coef(3,:));
    D = funD(coef(4,:));
end
%--------------------------Calculate BRL current-------------------------%

i = ((A*eta.^B).^-4+(C*eta.^D).^-4).^-0.25; % normalized current equation(8)

iIth = n*Jr*i;                              % BRL theoretical ion current 

%------------------------------------------------------------------------%

    
end

