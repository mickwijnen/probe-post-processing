%------------------------------------------------------------------------%
% Langmuir Probe Post-processing
% Comparison of LP and RFC-LP measurements
%
% Version:  1.0
% Data:     09/10/2016
% Author:   M. Wijnen
%------------------------------------------------------------------------%

function RFComparator
clc
close all

%% Input Parameters

RFpower = 100;        %  RF antenna input power

%-----------------------------File name----------------------------------%

filenameLP = 'KEITHLEY_20161005_105846_Z0X0Y0'; % LP filename
filenameRFCLP = 'KEITHLEY_20161005_110007_Z0X0Y0'; % RFCLP filename

probes = {filenameLP filenameRFCLP};

SAVE = 1;

%-------------------------File management--------------------------------%

date = datestr(now,'yyyymmdd');     % generates today's date string

savedir = fullfile('./results',date,'RF-compensation',...
    strcat(filenameLP(1:17),'_',num2str(RFpower),'W')); % save directory

if ~exist(savedir)
    mkdir(savedir);
end

warning('off','MATLAB:Axes:NegativeDataInLogAxis') % ignore logplot warning

%----------------------------File Import---------------------------------%

file =  char(fullfile('./raw data',strcat(probes(1),'.txt'))); 

inputdata = importFile(file); % imports data file to table

V   = inputdata.V;           % initialize voltage array

Iraw = inputdata.Iraw;        % initialize raw current array

SAT = find(V == 9999.999,1);  % remove saturated data points

if ~isempty(SAT)
    V = V(1:SAT-1);
    Iraw = Iraw(1:SAT-1);
end

file =  char(fullfile('./raw data',strcat(probes(2),'.txt'))); 

inputdata = importFile(file); % imports data file to table

RFV    = inputdata.V;           % initialize voltage array

RFIraw = inputdata.Iraw;        % initialize raw current array

SAT = find(RFV == 9999.999,1);  % remove saturated data points

if ~isempty(SAT)
    RFV = RFV(1:SAT-1);
    RFIraw = RFIraw(1:SAT-1);
end


%-------------------------------Smoothing--------------------------------%
    
span = 0.01;                            % define span for rloess algorithm 

I = smooth(V,Iraw,span,'rloess');% smooth data using rloess

sErr = 100*abs(I-Iraw)./abs(Iraw);   % smoothing error

RFI = smooth(RFV,RFIraw,span,'rloess');% smooth data using rloess

RFsErr = 100*abs(RFI-RFIraw)./abs(RFIraw);   % smoothing error

%--------------------------------Plotting--------------------------------%

IVplot = figure('Name', 'IV-curve')
plot(V,I)
hold on
grid on
plot(RFV,RFI)
title(['LP vs. RFCLP - IV curve @ ' num2str(RFpower) ' W'])
xlabel('Voltage [V]')
ylabel('Current [A]') 
vline(0,'k:')
hline(0,'k:')

legend('LP','RFCLP','Location','NorthWest')

Ionplot = figure('Name', 'Ion-current')
plot(V,I,'o','MarkerSize',3)
hold on
grid on
plot(RFV,RFI,'o','MarkerSize',3)
legend('LP','RFCLP','Location','SouthEast')
set(gca,'YLim',[-0.00005 0])
title(['LP vs. RFCLP - ion collection @ ' num2str(RFpower) ' W'])
xlabel('Voltage [V]')
ylabel('Current [A]') 

LogIeplot = figure('Name', 'Log Ie')
semilogy(V,I,'o','MarkerSize',3)
hold on
grid on
semilogy(RFV,RFI,'o','MarkerSize',3)
xlim([20 100])
legend('LP','RFCLP','Location','SouthEast')
title(['LP vs. RFCLP - electron collection @ ' num2str(RFpower) ' W'])
xlabel('Voltage [V]')
ylabel('Log of current') 

plotH = [IVplot Ionplot LogIeplot];
nameExt = {'_IV' '_ion' '_logIe'};

if SAVE
    for k = 1:length(plotH)
       saveas(plotH(k),...
         char(fullfile(savedir,strcat(filenameLP(1:17),'_',num2str(RFpower),'W', nameExt(k)))),'fig');
       saveas(plotH(k),...
         char(fullfile(savedir,strcat(filenameLP(1:17),'_',num2str(RFpower),'W', nameExt(k)))),'png');
    end
end


