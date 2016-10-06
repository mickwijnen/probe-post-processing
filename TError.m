function [Err] = TError(Var, V, I, iBounds, eBounds, const, probe)
% ERROR Generates the theoretical ion current and calculates the difference
%       with the measured current in the ion region. It then substracts the 
%       theoretical ion current from the measured current and compares this 
%       to the theoretical electron current in the electron region.
%       The sums of the difference vector are iErr and eErr. The total
%       error Err = iErr + eErr
%
%
% Output:           Err: sum of ion and electron fit error
%                   iIth: theoretical ion current
% Input:
%                   Var: all four variables, n iVs Te eVs
%                   V: measured voltage array
%                   I: measured current array
%                   iBounds: 1x2 vector with io fit bounds
%                   eBounds: 1x2 vector with electron fit bounds
%                   const: struct with natural constants e eps0 mi me 
%                   probe: struct with probe dimensions, Rp Lp S AR

iVar = Var(1:2);

eVar = Var(3:4);

[iErr, iIth] = iError(iVar, eVar, V, I, iBounds, const, probe);

[eErr, ~] = eError(eVar, iVar, V, I, iIth, eBounds, const, probe);

Err = iErr + eErr;

end


