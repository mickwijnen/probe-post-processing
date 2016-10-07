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

%-------------------------File management--------------------------------%

file =  fullfile('./raw data',strcat(filename,'.txt')); 

date = datestr(now,'yyyymmdd');     % generates today's date string

savedir = fullfile('./results',date,'Bohm',filename);  % specifies save directory

if ~exist(savedir)
    mkdir(savedir);
end

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
plot(V, Iraw,'o',V,I,'-','MarkerSize',4)       % plot raw vs smoothed data

ylabel('probe current [A]')

line(xlim, [0 0],'color','k','linewidth',0.5,'linestyle','--') %x-axis
line([0 0], ylim,'color','k','linewidth',0.5,'linestyle','--') %y-axis

grid minor
title('Comparison of raw and smoothed IV-curve')
legend('raw data', 'smoothed data','Location','NorthWest')

subplot(2,1,2)                         % plot relative smoothing error
plot(V, sErr,'o', V, ones(length(V))*mean(sErr),'r-','MarkerSize',3)

title('Smoothing error')
legend('error', sprintf('mean=%4.2f %%',mean(sErr)),'Location','NorthWest')

xlabel('Voltage (V)')
ylabel('diff/data [%]')
ylim([0 10])

%--------------------------Floating Potential----------------------------%

idxV0 = find(V>=0,1);           % find the index of V = 0
        
idxVf = find(I>0,1);            % find the index of the first I > 0

                                % interpolate the floating potential
Vf = V(idxVf-1) - I(idxVf-1)*(V(idxVf)-V(idxVf-1))/(I(idxVf)-I(idxVf-1));

dVf = 0.01;                  % relative error in Vf

%-----------------------------Derivatives--------------------------------%

I = -1000*I;                    % scale data positive and to mA range

dIdV = gradient(I,V);           % generate the first derivative dI/dV

%--------------------------Estimate Te and Vs----------------------------%

derivPlot = figure('Name','derivative','NumberTitle','off');
plot(V,dIdV);                     % plot derivative dI/dV

title('Estimation of V_s','Interpreter','tex')
xlabel('Voltage [V]')
ylabel('$\frac{dI}{dV}$','Interpreter','latex','Rotation',0,'FontSize',16)
      
set(gca,'XLim',[-20, max(V)])     % rescale window -20 to Vmax

display('Choose estimation:  Vs') % prompt user 

[x,~] = ginput(1);           % select x values of Vs

idxVs = knnsearch(V,x);      % find index of Vs
Vs = V(idxVs);               % set guess of Vs

vline(Vs,'k:')               % denote Vs on plot

dVs = 0.5;                   % error in Vs

annotation(derivPlot,'textbox',[0.5 0.1 0.25 0.1],...
    'String',['V_s = ' num2str(Vs,'%.1f') char(177) num2str(dVs) ' V'],...
    'FitBoxToText','on') 

%-----------------------------Ion Fitting--------------------------------%

I = -I/1000;                            % rescale I to [A]

ionPlot = figure('Name','Ion current','NumberTitle','off');

plot(V(1:idxV0),I(1:idxV0),'o','MarkerSize',3)   % plot ion current
hold on

xlim([V(1) 0])
xlabel('Voltage [V]')
ylabel('Current [A]')
title('Ion current')

display('Choose ion fit interval.')     % prompt user for fit interval
[x, ~] = ginput(2);

vline([min(x) max(x)],'k:');                    % denote fit boundaries

idxL = knnsearch(V,min(x));             % index of lower bound
idxU = knnsearch(V,max(x));             % index of upper bound

[iC,iS]= polyfit(V(idxL:idxU),I(idxL:idxU),1);% linear fit of ion fit range

iFitErr = sqrt(diag(inv(iS.R)*(iS.R')).*iS.normr.^2/iS.df);

Ii = iC(1)*V+iC(2);                       % extrapolate ion current

plot(V,Ii);                             % plot ion fit

legend('ion current','ion current fit','Location','SouthEast')

iDiff = abs(Ii-I)./abs(I);              % calculate relative difference

iErr  = mean(iDiff(idxL:idxU));         % calculate mean ion fit error       

%----------------------------Electron Fitting----------------------------%

Ie = I-Ii;                              % subtract ion current from total

electronPlot = figure('Name','Ion current','NumberTitle','off');

semilogy(V,Ie,'o',V,I,'MarkerSize',3)   % logplot electron current
hold on

vline(Vs,'k--')                         % denote plasma potential

xlim([-20 20])
xlabel('Voltage [V]')
ylabel('Log(I_e)','Interpreter','tex')
title('Log of electron current')

display('Choose electron fit interval.')% prompt user for fit interval
[x, ~] = ginput(2);

vline(x,'k:');                          % denote fit boundaries

idxL = knnsearch(V,min(x));             % lower fit boundary index
idxU = knnsearch(V,max(x));             % upper fit boundary index

[eC,eS] = polyfit(V(idxL:idxU),log(Ie(idxL:idxU)),1);% linear fit of log Ie

eFitErr = sqrt(diag(inv(eS.R)*(eS.R')).*eS.normr.^2/eS.df); % error in Te

Ief= exp(eC(1)*V + eC(2));                % electron current fit

semilogy(V,Ief)                           % plot fitted electron current

legend('electron current', 'total current', 'electron current fit',...
    'Location','SouthEast')

eDiff = abs(log(Ief)-log(I))./abs(log(I)); % calculate relative difference

eErr  = mean(eDiff(idxL:idxU));        % calculate mean electron fit error

%----------------------------Plasma Parameters---------------------------%

Te     = 1/eC(1);                     % calculate electron temperature
ReldTe = (eFitErr(1)/eC(1));          % relative error in Te
dTe    = ReldTe*Te;                      % absolute error in Te
    
annotation(electronPlot,'textbox',[0.15 0.8 0.25 0.1],...
    'String',['T_e = ' num2str(Te,'%.2f') char(177)...
    num2str(dTe,'%.2f') ' eV'],'FitBoxToText','on') 

Cs = sqrt(e*Te/mi);             % calculate Bohm velocity
    
Iis = iC(1)*Vf + iC(2);         % ion saturation = ion current at Vf 

ReldIis = sqrt((Vf*iC(1))^2*((iFitErr(1)/iC(1))^2+... % error in iIsat
    (dVf/Vf)^2)+iFitErr(2)^2);

ReldS = 0.1;                            % relative error in probe surface

n = -Iis/e/S/Cs/exp(-0.5);              % calculate the plasma density    

Reldn = sqrt(ReldIis^2 + (0.5*dTe)^2 + ReldS^2); % relative error in density

dn = n*Reldn;                            % absolute error in density

annotation(ionPlot,'textbox',[0.4 0.8 0.25 0.1],...
    'String',['n = (' num2str(n/1e17,'%.2f') char(177)...
    num2str(dn/1e17,'%.2f') ')\cdot10^{17} m^{-3}'],'FitBoxToText','on') 

DL = 1000*sqrt(eps0*Te/e/n);  % calculate the Debye length

xi = Rp/DL;                   % calculate the probe parameter

%--------------------------------Results---------------------------------%

results = fullfile(savedir,strcat(filename,'_results','.txt'));

fID = fopen(results,'w');

fprintf(fID,['Results for thin sheath approximation \r\n']);
fprintf(fID,['Plasma density = (' num2str(n/1e17,'%.2f') char(177)...
    num2str(dn/1e17,'%.2f') ')e+17 m^-3 \r\n']);
fprintf(fID,['Electron temperature = ' num2str(Te,'%.2f') char(177)...
    num2str(dTe,'%.2f') ' eV \r\n']);
fprintf(fID,['Plasma potential = ' num2str(Vs,'%.2f') char(177) num2str(dVs) ' V \r\n']);
fprintf(fID,['Floating potential = ' num2str(Vf,'%.1f') char(177) num2str(dVf) ' V \r\n']);
fprintf(fID,['Debye length = ' num2str(DL,'%.2f') ' mm \r\n']);
fprintf(fID,['Probe parameter = ' num2str(xi,'%.2f') '\r\n \r\n']);

fprintf(fID,['Ion fit error = ' num2str(iErr*100,'%.2f') ' %% \r\n']);
fprintf(fID,['Electron fit error = ' num2str(eErr*100,'%.2f') ' %% \r\n']);
fprintf(fID,['Total fit error = ' num2str(iErr+eErr,'%.2f') ' %% \r\n  \r\n']);

fclose(fID);

nameFig = [derivPlot ionPlot electronPlot];  % figure handles
nameExt = {'_dIdV' '_iFit' '_eFit'};       % figure name extensions

for k = 1:length(nameFig)               % save all figures as FIG and PNG
    saveas(nameFig(k),...
        char(fullfile(savedir,strcat(filename,nameExt(k)))),'fig');
    saveas(nameFig(k),...
        char(fullfile(savedir,strcat(filename,nameExt(k)))),'png');
end
end





