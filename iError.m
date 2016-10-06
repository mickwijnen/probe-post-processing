function [iErr, iIth] = iError(iVar, eVar, V, I, iBounds, const, probe)
%iError Generates a theoretical ion current and calculates the difference
%       with the measured current. Error is the sum of the absolute
%       difference vector.
%
% Output:           iErr: error in ion fit
%                   iIth: theoretical ion current
% Input:
%                   
%                   iVar: ion variables, n iVs
%                   eVar: electron variables, Te eVs
%                   V: measured voltage array
%                   I: measured current array
%                   iBounds: 1x2 vector with lower/upper index of interval
%                   const: struct with natural constants e eps0 mi me 
%                   probe: struct with probe dimensions, Rp Lp S AR

[iIth, DL, xi]  = paraBRL(V, iVar, eVar, const, probe); % theoretical curve

I = -I;                    % negate current data

% diff = abs(iIth-I);        % calculate difference between data and theory

% iErr = sum(diff(iBounds(1):iBounds(2))); % calculate error over interval 

reldiff = abs(iIth-I)./abs(I);

iErr = mean(reldiff(iBounds(1):iBounds(2)));

end

