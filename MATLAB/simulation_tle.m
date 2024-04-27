%% simulation using TLE's

clear, close all
%% Figure settings
set(0,'DefaultFigureColor',[1 1 1]); set(0,'DefaultLineLineWidth',0.9);
set(0,'DefaultAxesFontSize',12); set(0,'DefaultTextInterpreter','latex');
set(0,'DefaultAxesFontWeight','Normal');
set(0,'DefaultAxesTitleFontWeight','Bold');
%% load in info

load('data\obs_split.mat')

%% ODE settings
options = odeset('RelTol',1e-13,'AbsTol',1e-15);

%% Select Pass
tab = tabp{300}; % looking at pass 300 - Fylingdales SW 

%% Select Satellite
tab = tab(tab.fin_tag == 47492,:); % change as needed

%% get the  day frac to set up SGP4
% SGP-4 propagator is used to propagate TLEs. TLEs need to have fractional
% days
Tfrac = tab.day(1) + tab.hr(1)/24 + tab.min(1)/1440 + tab.s(1)/86400;
year = tab.year(1); 
% corresponds to 1st avail tle 
%% set up SGP4
tle_file = 'TLEs\47492_12.tle';
[sc, sat, sat_ID] = gnc.generate_sgp4_scenario(tab.time(1),tab.time(end),...
    1,tle_file,'sgp4');
[pos,vel] = states(sat,tab.time(1),'CoordinateFrame','inertial');
xm1 = [pos*1e-3;vel*1e-3];

%% initialize vars
% load('RadarInfo.mat') % sent

N = size(tab,1); % number of observations
x_sgp4 = nan(6,N);
x_sgp4(:,1) = xm1;
xode = nan(6,N);
xode(:,1) = x_sgp4(:,1);
Time = zeros(1,N);
y_real = nan(4,N);
y_exp = nan(4,N);
y_real(:,1) = [tab.range(1);tab.range_rate(1);tab.el(1);tab.az(1)];
snsr_idx = find(radar_table.sensor_num==tab.sensor_num(1)); % which radar is this meas from
lla = [radar_table.Latitude(snsr_idx),radar_table.Longitude(snsr_idx),...
    radar_table.Altitude(snsr_idx)]; % latitude, longitude, and altitude of the radar 
y_exp(:,1) = gnc.gen_observation_fn(xode(:,1),lla,tab.time(1));

% USE ODE45 INSTEAD
for i = 2:N
    %% sgp4 update
    [pos,vel] = states(sat,tab.time(i),'CoordinateFrame','inertial');
    x_sgp4(:,i) = [pos*1e-3;vel*1e-3];
end











