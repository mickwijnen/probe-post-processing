function dTe = TeUnCert(Var, V, eIth, Ie, eBounds)
%TEUNCERT Calculates the uncertainty in the electron temperature
%   The electron temperature is found by minimizing abs(log(Ie)-log(eIth)
%   which is the equivalent of using a least squares linear fit. 
%   The error in the first fit constant C1 of a LSQ is given by the square
%   root of SSxy2/SSxx (defined below). Since Te = 1/C1 this can be used to
%   calculate dTe the error in Te.
%
%   Output: dTe: error in Te
%
%   Input: 
%           V: array of voltage data
%           eIth: theoretical electron current
%           Ie: measured electron current
%           Te: electron temperature
%           eBounds: indices of electron fit interval

Te = Var(3);

V    = V(eBounds(1):eBounds(2));        % subset of V over fit interval
Ie   = Ie(eBounds(1):eBounds(2));       % subset of Ie over fit interval
eIth = eIth(eBounds(1):eBounds(2));     % subset of eIth over fit interval

SSxy2 = sum((log(Ie)-log(eIth)).^2)/(length(Ie)-2); % variance of log of Ie and eIth
SSxx  = sum((V-mean(V)).^2);                % mean of V

Del = sqrt(SSxy2/SSxx);                    % relative error in hypothetical
                                           % linear fit constant
dTe = Del*Te^2;                            % absolute error in Te

end

