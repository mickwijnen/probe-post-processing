function [eErr, eIth, Ie] = eError(eVar, iVar, V, I, iIth, eBounds, const, probe)
%eError Generates theoretical electron current and calculates the error
%       with the measured current. Error is the sum of the absolute
%       difference vector of the natural log of te current.
%
% Takes as input:
%                   V: measured voltage array
%                   I: measured current array
%                   iVar: ion variables, n iVs
%                   eVar: electron variables, Te eVs
%                   iIth: theoretical ion current
%                   eBounds: 1x2 vector with lower/upper index of interval
%                   const: struct with natural constants e eps0 mi me 
%                   probe: struct with probe dimensions, Rp Lp S AR

e  = const.e;       % [C] elementary charge
me = const.me;      % [kg] electron mass
S  = probe.S;       % [mm^2] probe surface area

n   = iVar(1)*1e17; % [m^-3] plasma density
Te  = eVar(1);      % [eV] electron temperature
eVs = eVar(2);      % [V] eleron plasma potential

Ce = sqrt(e*Te/2/pi/me);            % [m/s] electron thermal current

eIth = e*n*S*Ce*exp(-(eVs-V)/Te);    % theoretical electron current
    
Ie = I + iIth;                      % measured electron current corrected
                                    % for theoretical ion current
                                   
Ie = Ie.*(Ie > 0);                  % set all negative values zero
                                  
diff = abs(log(Ie)-log(eIth));      % error of electron fit

eErr = sum(diff(eBounds(1):eBounds(2))); % limit error to eBounds

end

