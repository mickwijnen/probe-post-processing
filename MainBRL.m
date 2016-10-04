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
clear all

%% Input Parameters
%--------------------------Probe properties------------------------------%

Rp =  0.00635;  %[mm] probe radius 
Lp =  5;        %[mm] probe length 

%---------------------------Gas properties-------------------------------%

AN = 39.95;     %[amu] atomic mass number 

%-----------------------------File name----------------------------------%

filename = 'KEITHLEY_20160930_153846_Z0X0Y0'; %filename WITHOUT extension

%% Preprocessing
%-----------------------------Constants----------------------------------%

e   = 1.60217662e-19;    %[C]   electron charge
me  = 9.109382910e-31;   %[kg]  electron mass
ep0 = 8.854187817e-12;   %[F/m] Vacuum permittivity
amu = 1.660539040e-27;   %[kg]  atomic mass unit

S   = 2*pi*(Rp*1e-3)*(Lp*1e-3); %[m^2]  probe surface
AR  = 1+Rp/2/Lp;                %[]     aspect ratio
mi  = AN*amu;                   %[kg]   ion mass

% structure
const = struct('e', e,'me', me,'mi', mi, 'ep0', ep0, 'S', S, 'AR', AR); 

%-------------------------File management--------------------------------%

file =  strcat(filename,'.txt');

date = datestr(now,'yyyymmdd');

savedir  =  fullfile('./results',date,filename)












end