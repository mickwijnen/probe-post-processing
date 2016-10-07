function dn = nUnCert(Var, iErr, dTe, probe)
%nUnCert calculates the error in the plasma density
%   Detailed explanation goes here

n = Var(1);

Te = Var(3);

ReldS  = 0.1;           % relative error in probe surface 

ReldTe = (dTe/Te);      % relative error in electron temperature

Reldi  = 0.03;          % relative error in Chen's parametrization

%calculate absolute error
dn  =  n*sqrt(ReldS^2 + ReldTe^2 + Reldi^2 + iErr^2); 

end

