%------------------------------------------------------------------------%
% Langmuir Probe Post-processing
% Based on F.F. Chen's iterative adapatation of BRL theory
%
% Version:  1.0
% Data:     04/10/2016
% Author:   M. Wijnen
%------------------------------------------------------------------------%

function MainBRL
clc
close all

%% Input Parameters
%--------------------------Probe properties------------------------------%

Rp =  0.127;  %[mm] probe radius 
Lp =  5;        %[mm] probe length 

%---------------------------Gas properties-------------------------------%

AN = 39.95;     %[amu] atomic mass number 

%-----------------------------File name----------------------------------%

filename = 'KEITHLEY_20160930_153846_Z0X0Y0'; %filename WITHOUT extension

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

file =  strcat(filename,'.txt');    % appends file extension to filename

date = datestr(now,'yyyymmdd');     % generates today's date string

savedir = fullfile('./results',date,filename);  % specifies save directory

%----------------------------File Import---------------------------------%

inputdata = importFile(file); % imports data file to table

V    = inputdata.V;           % initialize voltage array
Iraw = inputdata.Iraw;        % initialize raw current array

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

IdVdI = I./dIdV;                % generate the quotient of I to dI/dV

IdVdI = (IdVdI>0).*IdVdI;

%--------------------------Estimate Te and Vs----------------------------%

derivfig = figure('Name','derivative','NumberTitle','off');
[hAx,hdIdV,hIdVdI] = plotyy(V,dIdV,V,IdVdI);

title('Estimation of T_e and V_s','Interpreter','tex')
xlabel('Voltage [V]')
set(hAx,'XLim',[0, 15])

ylabel(hAx(1),'$\frac{dI}{dV}$', ...
'Interpreter','latex','Rotation',0,'FontSize',16)
ylabel(hAx(2),'$T_e [eV]$','Interpreter','latex','FontSize',12)

display('Choose estimation:  Vs(1st) and  Te(2nd)') % prompt user 

[x,~] = ginput(2);              % select x values of Vs and Te

idxVs = knnsearch(V,x(1));      % find index of Vs
gVs = V(idxVs);                 % set guess of Vs

iVs = gVs;                      % set iVs to guess of Vs
eVs = gVs;                      % set eVs to guess of Vs

Te = IdVdI(knnsearch(V,x(2)));  % set guess of Te

%----------------------------Estimate of n-------------------------------%
I = -I/1000;                          % rescale I to [A]

idxU = knnsearch(V,V(1)+50);            % upper idx of fit range: Vmin+50V

C = polyfit(V(1:idxU),I(1:idxU),1);     % linear fit of ion fit range

Cs = sqrt(e*Te/mi);                     % Bohm velocity

n = (C(1)*Vf+C(2))/(-0.61*e*S*Cs); % estimate density with Bohm

iVar = [n/1e17 iVs];
eVar = [Te eVs];

%-------------------------Default fit interval---------------------------%

iVmin = V(V == V(1)+10);    % default ion fit lower limit
iVmax = 0;                  % default ion fit upper limit

eVmin = eVs-5*Te;           % default electron fit lower limit
eVmax = eVs-Te/2;           % default electron fit upper limit

iLidx = knnsearch(V,iVmin); % index of ion fit lower limit
iUidx = knnsearch(V,iVmax); % index of ion fit upper limit

iBounds  = [iLidx iUidx];
iVbounds = [iVmin iVmax];

eLidx = knnsearch(V,eVmin); % index of electron fit upper limit
eUidx = knnsearch(V,eVmax); % index of electron fit upper limit

eBounds =  [eLidx eUidx];
eVbounds = [eVmin eVmax];

%-----------------------Plot Results from Estimate-----------------------%

% create plot handle
Isquare = figure('Name','I squared','NumberTitle','off');

% generate theoretical ion curve with new iVar values
[iIth, ~] = paraBRL(V, iVar, eVar, const, probe);

% plot the theoretical and measured squared ion currents
Isqrplot(Isquare, V, I, iVar, eVar, iVbounds, const, probe);

% create plot handle
TeFit = figure('Name','Electron Temperature','NumberTitle','off');

%----------------------------Optimize Fits-------------------------------%

% set plot options
options = optimset('Display','none','PlotFcns',@optimplotx);

for i = 1:2
% run optimizer with iError to minimize iErr
[iVar, iErr] = fminsearch(@(iVar)...
iError(iVar, eVar, V, I, iBounds, const, probe), [n/1e17 iVs], options);

% generate theoretical ion curve with new iVar values
[iIth, ~] = paraBRL(V, iVar, eVar, const, probe);

% plot the theoretical and measured squared ion currents
Isqrplot(Isquare, V,I, iVar, eVar, iVbounds, const, probe);

if i == 1
[iBounds, iVbounds] = iSetfitbnd(V); % prompt user to define new fit bounds
end

% run optimizer with eError to minimize eErr
[eVar, eErr] = fminsearch(@(eVar)...
eError(eVar, iVar, V, I, iIth, eBounds, const, probe), [Te eVs], options);

% Logplot of measured and theoretical electron currents
TeFitplot(TeFit,V,I,iIth,iVar,eVar,eVbounds,eBounds,const,probe)

if i == 1
[eBounds, eVbounds] = eSetfitbnd(V); % prompt user to define new fit bounds
end
end

iErr
eErr






end

%% Function for plotting squared currents

function [] = Isqrplot(Isquare, V,I, iVar, eVar, iVbounds, const, probe)

[iIth, ~] = paraBRL(V, iVar, eVar, const, probe);

Isqr = (I.*(I < 0)).^2;                    % square of current where I < 0

figure(Isquare)

plot(V,Isqr,'o',V,iIth.^2,'MarkerSize',4)   % plot data vs theory, squared
vline(iVbounds,'k:')                        % plot fit boundaries

title('Theoretical vs. measured square of ion current')
xlabel('Voltage [V]')
ylabel('Current squared [A^2]','Interpreter','tex')
legend('Measured current','Theoretical current')

end

%% Function for plotting electron fit

function [] = TeFitplot(TeFit,V,I,iIth,iVar,eVar,eVbounds, eBounds,const,probe)

figure(TeFit)

[eErr, eIth, Ie] = eError(eVar, iVar, V, I, iIth, eBounds, const, probe);

Vlow = V(find(Ie == 0,1,'last'));             % find lowest V where Ie > 0 

semilogy(V,Ie,'o',V,I,V,eIth,'MarkerSize',4); % logplot of electron current

limV = get(gca, 'XLim');                      % get current x limit
set(gca, 'Xlim',[Vlow limV(2)]);              % set lower x limit

vline(eVbounds,'k:')                          % plot fit boundaries

title('Natural log of theoretical vs. measured electron current')
xlabel('Voltage [V]')
ylabel('Log(I_e)','Interpreter','tex')
legend('electron current','total current','theoretical electron current',...
'Location','SouthEast')

end

%% function to set ion fit boundaries

function [iBounds, iVbounds] = iSetfitbnd(V)

    display('Choose new ion fit boundaries.')
    [x,~] = ginput(2);
    
    if length(x) == 2;
    
        iVmin = min(x);
        iVmax = max(x);
        display(strcat({'New ion fit boundaries are '},...
        {num2str(iVmin)},{' and '}, {num2str(iVmax)},{' V'}))

    else
        display('Ion fit boundaries unchanged')
    end    
    
    iLidx = knnsearch(V,iVmin); % index of ion fit lower limit
    iUidx = knnsearch(V,iVmax); % index of ion fit upper limit

    iBounds  = [iLidx iUidx];
    iVbounds = [iVmin iVmax];
    
end

%% function to set electron fit boundaries
function [eBounds, eVbounds] = eSetfitbnd(V)

    display('Choose new electron fit boundaries.')
    [x,~] = ginput(2);

    if length(x) == 2;
    
        eVmin = min(x);
        eVmax = max(x);
        display(strcat({'New electron fit boundaries are '},...
        {num2str(eVmin)},{' and '}, {num2str(eVmax)},{' V'}))

    else
        display('Electron fit boundaries unchanged')
    end
    
    eLidx = knnsearch(V,eVmin); % index of electron fit upper limit
    eUidx = knnsearch(V,eVmax); % index of electron fit upper limit

    eBounds =  [eLidx eUidx];
    eVbounds = [eVmin eVmax];
    
end

%% function to generate boundary parameters from (i/e)Vmin (i/e)Vmax
