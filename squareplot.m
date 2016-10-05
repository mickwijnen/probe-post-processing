function [] = Isqrplot(V,I,iIth)
%SQUAREPLOT Generate a plot of the measured vs. theoretical current squared
%   Detailed explanation goes here

Isqr = (I.*(I < 0)).^2;                    % square of current where I < 0

Isquare = figure('Name','I squared','NumberTitle','off');

plot(V,Isqr,'o',V,iIth.^2,'MarkerSize',4)  % plot data vs theory, squared

title('Theoretical vs. measured square of ion current')
xlabel('Voltage [V]')
ylabel('Current squared [A^2]','Interpreter','tex')

end

