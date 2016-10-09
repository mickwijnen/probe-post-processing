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
Lp =  5;      %[mm] probe length 

%---------------------------Gas properties-------------------------------%

AN = 39.95;     %[amu] atomic mass number 

%-----------------------------File name----------------------------------%

filename = 'KEITHLEY_20161005_113323_Z0X0Y0'; %filename WITHOUT extension

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

savedir = fullfile('./results',date,'BRL',filename);  % specifies save directory

if ~exist(savedir)
    mkdir(savedir);
end 

warning('off','MATLAB:Axes:NegativeDataInLogAxis') % ignore logplot warning

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
plot(V, Iraw,'o',V,I,'-','MarkerSize',3)

ylabel('probe current [A]')

line(xlim, [0 0],'color','k','linewidth',0.5,'linestyle','--') %x-axis
line([0 0], ylim,'color','k','linewidth',0.5,'linestyle','--') %y-axis

grid minor
title('Comparison of raw and smoothed IV-curve')
legend('raw data', 'smoothed data','Location','NorthWest')

subplot(2,1,2)
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
Vf = V(idxVf-1) - I(idxVf-1)*(V(idxVf)-V(idxVf-1))/(I(idxVf)-I(idxVf-1))

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

Vlow = get(gca, 'Xlim');
set(hAx,'XLim',[-20, Vlow(2)])
set(hAx(2),'Ylim',[0 8])
set(hAx(2),'YTick',[0:8])

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

idxU = knnsearch(V,V(1)+50);          % upper idx of fit range: Vmin+50V

C = polyfit(V(1:idxU),I(1:idxU),1);   % linear fit of ion fit range

Cs = sqrt(e*Te/mi);                   % Bohm velocity

n = (C(1)*Vf+C(2))/(-0.61*e*S*Cs);    % estimate density with Bohm

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
Isqr(Isquare, V, I, iVar, eVar, iVbounds, const, probe);

% create plot handle
LogTe = figure('Name','Electron Temperature','NumberTitle','off');

%----------------------------Optimize Fits-------------------------------%

% set plot options
options = optimset('Display','none','PlotFcns',@optimplotx);

for i = 1:3
% run optimizer with iError to minimize iErr
[iVar, iErr] = fminsearch(@(iVar)...
iError(iVar, eVar, V, I, iBounds, const, probe), [n/1e17 iVs], options);

% generate theoretical ion curve with new iVar values
[iIth, ~] = paraBRL(V, iVar, eVar, const, probe);

% plot the theoretical and measured squared ion currents
Isqr(Isquare, V,I, iVar, eVar, iVbounds, const, probe);

[iBounds, iVbounds] = iSetfitbnd(V); % prompt user to define new fit bounds

% run optimizer with eError to minimize eErr
[eVar, eErr] = fminsearch(@(eVar)...
eError(eVar, iVar, V, I, iIth, eBounds, const, probe), [Te eVs], options);

% Logplot of measured and theoretical electron currents
[Ie eIth] = TeFit(LogTe,V,I,iIth,iVar,eVar,eVbounds,eBounds,const,probe);

[eBounds, eVbounds] = eSetfitbnd(V); % prompt user to define new fit bounds

end

Var = [iVar eVar]           % combine all variables

Err = iErr + eErr;          % total fit error

dTe = TeUnCert(Var, V, eIth, Ie, eBounds); % calculate Te error

dn = nUnCert(Var, iErr, dTe, probe); % calculate density error

annotation(Isquare,'textbox',[0.25 0.15 0.25 0.1],...
  'String',['n = (' num2str(Var(1),'%.2f') char(177) num2str(dn,'%.2f')...
  ')\cdot10^{17} m^{-3}'],'FitBoxToText','on')    % create annotation

annotation(LogTe,'textbox',[0.65 0.3 0.25 0.1],...
 'String',['Te = ' num2str(Var(3),'%.2f') char(177) num2str(dTe,'%.2f')...
 ' eV'],'FitBoxToText','on')                % create annotation
%----------------------------Final Optimization--------------------------%                 

[Var, OPTErr] = fminsearch(@(Var)...           % run final optimization
TError(Var, V, I, iBounds, eBounds, const, probe), Var, options);

[OPTErr, OPTiErr, OPTeErr] = TError(Var,V,I,iBounds,eBounds,const,probe);

OPTiVar = Var(1:2);            % set optimal ion variables
OPTeVar = Var(3:4);            % set optimal electron variables 

n  = Var(1)*1e17;
Te = Var(3);

DL = 1000*sqrt(eps0*Te/e/n); % [mm] Debye length
xi = Rp/DL;                  % Probe radius - Debye length ratio

%-------------------------------Results----------------------------------%

IsquareOPT = figure('Name','Optimized Density Fit','NumberTitle','off');

LogTeOPT = figure('Name','Optimized Electron Fit','NumberTitle','off');

[iIth] = Isqr(IsquareOPT,V,I,OPTiVar,OPTeVar,iVbounds,const,probe);

[Ie eIth] = TeFit(LogTeOPT,V,I,iIth,OPTiVar,...
    OPTeVar,eVbounds,eBounds,const,probe);

OPTdTe = TeUnCert(Var, V, eIth, Ie, eBounds); % calculate Te error

OPTdn = nUnCert(Var, OPTiErr, OPTdTe, probe); % calculate density error

annotation(IsquareOPT,'textbox',[0.25 0.15 0.25 0.1],...
'String',['n = (' num2str(Var(1),'%.2f') char(177) num2str(OPTdn,'%.2f')...
')\cdot10^{17} m^{-3}'],'FitBoxToText','on')    % create annotation

annotation(LogTeOPT,'textbox',[0.65 0.3 0.25 0.1],...
 'String',['Te = ' num2str(Var(3),'%.2f') char(177) num2str(OPTdTe,'%.2f')...
 ' eV'],'FitBoxToText','on')                % create annotation

% Create a text file to write results

results = fullfile(savedir,strcat(filename,'_results','.txt'));

fID = fopen(results,'w');

fprintf(fID,['Results for last consecutive optimization \r\n']);
fprintf(fID,['Plasma density = (' num2str(iVar(1), '%.2f') char(177)...
    num2str(dn,'%.2f') ')e+17 m^-3 \r\n']);
fprintf(fID,['Electron temperature = ' num2str(eVar(1),'%.2f') char(177)...
    num2str(dTe,'%.2f') ' eV \r\n']);
fprintf(fID,['Electron plasma potential = ' num2str(eVar(2),'%.2f') ' V \r\n']);
fprintf(fID,['Ion plasma potential = ' num2str(iVar(2),'%.2f') ' V \r\n \r\n']);

fprintf(fID,['Ion fit error = ' num2str(iErr*100,'%.2f') ' %% \r\n']);
fprintf(fID,['Electron fit error = ' num2str(eErr*100,'%.2f') ' %% \r\n']);
fprintf(fID,['Total fit error = ' num2str(Err*100,'%.2f') ' %% \r\n  \r\n']);

fprintf(fID,['Results for last consecutive optimization \r\n \r\n']);
fprintf(fID,['Plasma density = (' num2str(OPTiVar(1),'%.2f')...
    char(177) num2str(OPTdn,'%.2f') ')e+17 m^-3 \r\n']);
fprintf(fID,['Electron temperature = ' num2str(OPTeVar(1),'%.2f')... 
    char(177) num2str(OPTdTe,'%.2f') ' eV \r\n']);
fprintf(fID,['Electron plasma potential = ' num2str(OPTeVar(2),'%.2f') ' V \r\n']);
fprintf(fID,['Ion plasma potential = ' num2str(OPTiVar(2),'%.2f') ' V \r\n']);
fprintf(fID,['Floating potential = ' num2str(Vf,'%.2f') ' V \r\n']);
fprintf(fID,['Debye length = ' num2str(DL,'%.2f') ' mm \r\n']);
fprintf(fID,['Probe parameter = ' num2str(xi,'%.2f') '\r\n \r\n']);

fprintf(fID,['Ion fit error = ' num2str(OPTiErr*100,'%.2f') ' %% \r\n']);
fprintf(fID,['Electron fit error = ' num2str(OPTeErr*100,'%.2f') ' %% \r\n']);
fprintf(fID,['Total fit error = ' num2str(OPTErr*100,'%.2f') ' %% \r\n']);

fclose(fID);

nameFig = [Isquare, LogTe, IsquareOPT, LogTeOPT,];  % figure handles
nameExt = {'_Isqr' '_LogTe' '_IsqrOPT' '_LogTeOPT'};   % figure name extensions

for k = 1:length(nameFig)               % save all figures as FIG and PNG
    saveas(nameFig(k),...
        char(fullfile(savedir,strcat(filename,nameExt(k)))),'fig');
    saveas(nameFig(k),...
        char(fullfile(savedir,strcat(filename,nameExt(k)))),'png');
end


end

%% Function for plotting squared currents

function [iIth] = Isqr(Isquare,V,I,iVar,eVar,iVbounds,const,probe)

[iIth, ~] = paraBRL(V, iVar, eVar, const, probe);

Isqr = (I.*(I < 0)).^2;                    % square of current where I < 0

figure(Isquare)

plot(V,Isqr,'o',V,iIth.^2,'MarkerSize',3)   % plot data vs theory, squared
vline(iVbounds,'k:')                        % plot fit boundaries

title('Theoretical vs. measured square of ion current')
xlabel('Voltage [V]')
ylabel('Current squared [A^2]','Interpreter','tex')
legend('measured current','theoretical current')



end

%% Function for plotting electron fit

function [Ie, eIth] = TeFit(LogTe,V,I,iIth,iVar,eVar,...
          eVbounds,eBounds,const,probe)

figure(LogTe)

[eErr, eIth, Ie] = eError(eVar, iVar, V, I, iIth, eBounds, const, probe);

Vlow = V(find(Ie == 0,1,'last'));             % find lowest V where Ie > 0 
if isempty(Vlow)                              % if Ie > 0 everywhere find 
   Vlow = V(min(Ie))                          % find V where Ie = min(Ie)
end

semilogy(V,Ie,'o',V,I,V,eIth,'MarkerSize',3); % logplot of electron current

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

    persistent iVbnd iBnd;
    
    display('Choose new ion fit boundaries.')
    [x,~] = ginput(2);
    
    if length(x) == 2;
    
        iVmin = min(x);
        iVmax = max(x);
        display(strcat({'New ion fit boundaries are '},...
        {num2str(iVmin)},{' and '}, {num2str(iVmax)},{' V'}))
    
        iLidx = knnsearch(V,iVmin); % index of ion fit lower limit
        iUidx = knnsearch(V,iVmax); % index of ion fit upper limit

        iBounds  = [iLidx iUidx];
        iVbounds = [iVmin iVmax];
        
        iVbnd = iVbounds;
        iBnd = iBounds;
    
    else
        display('Ion fit boundaries unchanged')
        
        iVbounds = iVbnd;
        iBounds = iBnd;
    end    
    
  
end

%% function to set electron fit boundaries
function [eBounds, eVbounds] = eSetfitbnd(V)
 
    display('Choose new electron fit boundaries.')
    [x,~] = ginput(2);

    persistent eVbnd eBnd;
    
    if length(x) == 2;
    
        eVmin = min(x);
        eVmax = max(x);
        display(strcat({'New electron fit boundaries are '},...
        {num2str(eVmin)},{' and '}, {num2str(eVmax)},{' V'}))
    
        eLidx = knnsearch(V,eVmin); % index of electron fit upper limit
        eUidx = knnsearch(V,eVmax); % index of electron fit upper limit

        eBounds =  [eLidx eUidx];
        eVbounds = [eVmin eVmax];
        
        eVbnd = eVbounds;
        eBnd = eBounds;
    else
        display('Electron fit boundaries unchanged')
        
        eVbounds = eVbnd;
        eBounds = eBnd;
    end
          
end
%% function to calculate the electron temperature and plasma density errors

