function veloc = velocities(Te,const)
%VELOCITIES Calculates velocities based on Te
%   Calculates the ion sound speed Cs
%   Calculates the electron thermal velocity Ce
%   Calculates the ion thermal velocity Va

e = const.e;
mi = const.mi;
me = const.me;

Cs = sqrt(e*Te/mi);
Ce = sqrt(e*Te/2/pi/me);
Va = sqrt(e*Te/2/pi/mi);

veloc = struct('Cs', Cs, 'Ce', Ce, 'Va', Va);

end

