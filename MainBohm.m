%------------------------------------------------------------------------%
% Langmuir Probe Post-processing
% Based on thin sheath approximation and Bohm criterion
%
% Version:  1.0
% Data:     06/10/2016
% Author:   M. Wijnen
%------------------------------------------------------------------------%

function MainBohm
clc
close all

%% Input Parameters
%--------------------------Probe properties------------------------------%

Rp =  0.127;  %[mm] probe radius 
Lp =  5;        %[mm] probe length 

%---------------------------Gas properties-------------------------------%

AN = 39.95;     %[amu] atomic mass number 

%-----------------------------File name----------------------------------%

filename = 'KEITHLEY_20161005_112613_Z0X0Y0'; %filename WITHOUT extension

%% Pre-processing
%-----------------------------Constants----------------------------------%

e   = 1.60217662e-19;    %[C]   electron charge
me  = 9.109382910e-31;   %[kg]  electron mass
eps0 = 8.854187817e-12;   %[F/m] Vacuum permittivity
amu = 1.660539040e-27;   %[kg]  atomic mass unit

S   = 2*pi*(Rp*1e-3)*(Lp*1e-3); %[m^2]  probe surface
AR  = 1+Rp/2/Lp;                %[]     aspect ratio
mi  = AN*amu;                   %[kg]   ion mass

% structure
const = struct('e', e,'me', me,'mi', mi, 'eps0', eps0); 
probe = struct('Rp', Rp, 'Lp', Lp, 'S', S, 'AR', AR);

%-------------------------File management--------------------------------%

file =  fullfile('./raw data',strcat(filename,'.txt')); 

date = datestr(now,'yyyymmdd');     % generates today's date string

savedir = fullfile('./results',date,'Bohm',filename);  % specifies save directory

mkdir(savedir);

%----------------------------File Import---------------------------------%

inputdata = importFile(file); % imports data file to table

V    = inputdata.V;           % initialize voltage array

Iraw = inputdata.Iraw;        % initialize raw current array

SAT = find(V == 9999.999,1);  % remove saturated data points
if ~isempty(SAT)
    V = V(1:SAT-1);
    Iraw = Iraw(1:SAT-1);
end

%-----------------------------Smoothing----------------------------------%

span = 0.01;                        % define span for rloess algorithm 

I = smooth(V,Iraw,span,'rloess');   % smooth data using rloess

sErr = 100*abs(I-Iraw)./abs(Iraw);  % smoothing error
    
I = I/AR;                           % scale data with aspect ratio
    
%--------------------------Smoothed data plot----------------------------%
smoothfig = figure('Name','smoothed IV-curve','NumberTitle','off');

subplot(2,1,1)
plot(V, Iraw,'o',V,I,'-','MarkerSize',4)

ylabel('probe current [A]')

line(xlim, [0 0],'color','k','linewidth',0.5,'linestyle','--') %x-axis
line([0 0], ylim,'color','k','linewidth',0.5,'linestyle','--') %y-axis

grid minor
title('Comparison of raw and smoothed IV-curve')
legend('raw data', 'smoothed data','Location','NorthWest')

subplot(2,1,2)
plot(V, sErr,'o', V, ones(length(V))*mean(sErr),'r-','MarkerSize',4)

title('Smoothing error')
legend('error', sprintf('mean=%4.2f %%',mean(sErr)),'Location','NorthWest')

xlabel('Voltage (V)')
ylabel('diff/data [%]')
ylim([0 10])

%--------------------------Floating Potential----------------------------%

idxV0 = find(V>=0,1);           % find the index of V = 0
        
idxVf = find(I>0,1);            % find the index of the floating potential

if I(idxVf-1)<0                 % interpolate the floating potential
    Vf = (V(idxVf)+ V(idxVf-1))/2;   
else
    Vf = (V(idxVf)+ V(idxVf+1))/2;
end 

%-----------------------------Derivatives--------------------------------%

I = -1000*I;                    % scale data positive and to mA range

dIdV = gradient(I,V);           % generate the first derivative dI/dV

%--------------------------Estimate Te and Vs----------------------------%

derivfig = figure('Name','derivative','NumberTitle','off');
hAx = plot(V,dIdV);

title('Estimation of V_s','Interpreter','tex')
xlabel('Voltage [V]')
ylabel('$\frac{dI}{dV}$','Interpreter','latex','Rotation',0,'FontSize',16)

Vlow = get(gca, 'Xlim');
set(gca,'XLim',[-20, Vlow(2)])

display('Choose estimation:  Vs') % prompt user 

[x,~] = ginput(1);              % select x values of Vs and Te

idxVs = knnsearch(V,x);      % find index of Vs
Vs = V(idxVs);               % set guess of Vs

%----------------------------Estimate of n-------------------------------%

I = -I/1000;                            % rescale I to [A]

ionPlot = figure('Name','Ion current','NumberTitle','off');

plot(V(1:idxV0),I(1:idxV0),'o','MarkerSize',3)             % plot ion current
hold on

xlim([V(1) 0])
xlabel('Voltage [V]')
ylabel('Current [A]')
title('Ion current')

display('Choose ion fit interval.')     % prompt user to choose fit interval
[x, ~] = ginput(2)

Vl = min(x);                            % set lower fit boundary
Vu = max(x);                            % set upper fit boundary

idxL = knnsearch(V,Vl);                 % index of lower bound
idxU = knnsearch(V,Vu);                 % index of upper bound

C = polyfit(V(idxL:idxU),I(idxL:idxU),1)% linear fit of ion fit range

Ii = C(1)*V+C(2);                       % extrapolate ion current

plot(V,Ii);       % plot ion fit

Ie = I-Ii;                              % subtract ion current from total

electronPlot = figure('Name','Ion current','NumberTitle','off');

semilogy(V,Ie,'o','MarkerSize',3)
hold on

xlim([-20 20])
xlabel('Voltage [V]')
ylabel('Log(I_e)','Interpreter','tex')
title('Log of electron current')

display('Choose electron fit interval.')     % prompt user to choose fit interval
[x, ~] = ginput(2)

idxL = knnsearch(V,min(x));
idxU = knnsearch(V,max(x));

Q = polyfit(V(idxL:idxU),log(Ie(idxL:idxU)),1)

Ies = exp(Q(1)*V + Q(2));

semilogy(V,Ies)

Te = 1/Q(1)

Cs = sqrt(e*Te/mi);

Iis = C(1)*Vf + C(2)

n = -Iis/e/S/Cs/0.61

DL = 1000*sqrt(eps0*Te/e/n)

xi = Rp/DL



