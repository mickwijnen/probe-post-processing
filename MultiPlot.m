%------------------------------------------------------------------------%
% Langmuir Probe Post-processing
% Comparison of LP and RFC-LP measurements
%
% Version:  1.0
% Data:     09/10/2016
% Author:   M. Wijnen
%------------------------------------------------------------------------%

function MultiPlot
clc
close all

%% Input Parameters

%-----------------------------File name----------------------------------%

files    = {
'KEITHLEY_20161005_105846_Z0X0Y0' 
'KEITHLEY_20161005_110243_Z0X0Y0' 
'KEITHLEY_20161005_110831_Z0X0Y0' 
'KEITHLEY_20161005_111337_Z0X0Y0' 
'KEITHLEY_20161005_111511_Z0X0Y0'
'KEITHLEY_20161005_111842_Z0X0Y0'
'KEITHLEY_20161005_112519_Z0X0Y0'} 

%-------------------------File management--------------------------------%

date = datestr(now,'yyyymmdd');     % generates today's date string

filename = char(files(1));
filename = filename(1:17);

savedir = fullfile('./results', date, filename); % save directory

if ~exist(savedir)
    mkdir(savedir);
end

warning('off','MATLAB:Axes:NegativeDataInLogAxis') % ignore logplot warning

%----------------------------=-------------------------------------------%

IVcurve = figure('Name','IVcurve') 

for i = 1:length(files)
    
file =  char(fullfile('./raw data',strcat(files(i),'.txt'))); 

inputdata = importFile(file); % imports data file to table

V   = inputdata.V;           % initialize voltage array

Iraw = inputdata.Iraw;        % initialize raw current array

SAT = find(V == 9999.999,1);  % remove saturated data points

if ~isempty(SAT)
    V = V(1:SAT-1);
    Iraw = Iraw(1:SAT-1);
end

span = 0.01;                        % define span for rloess algorithm 

I = smooth(V,Iraw,span,'rloess');   % smooth data using rloess

sErr = 100*abs(I-Iraw)./abs(Iraw);  % smoothing error

plot(V,I)
hold on

legendInfo{i} = [num2str(50+50*i) 'W']; % or whatever is appropriate

end

legend(legendInfo,'Location','NorthEast');
