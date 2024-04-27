%% Confirm ODE45 Dynamics
% spacex depluyment

%% Clearing
clear, close all
%% Figure settings
set(0,'DefaultFigureColor',[1 1 1]); set(0,'DefaultLineLineWidth',0.9);
set(0,'DefaultAxesFontSize',12); set(0,'DefaultTextInterpreter','latex');
set(0,'DefaultAxesFontWeight','Normal');
set(0,'DefaultAxesTitleFontWeight','Bold');
%% load in info
dir_name = 'C:\Users\lpdav\Documents\CU\ODR\SpaceXData';
load([dir_name,'\obs_split'])

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
addpath('C:\Users\lpdav\Documents\CU\ODR\SpaceXData\tles')
tle_file = '47492_12.tle';
[sc, sat, sat_ID] = generate_sgp4_scenario(tab.time(1),tab.time(end),...
    1,tle_file,'sgp4');
[pos,vel] = states(sat,tab.time(1),'CoordinateFrame','inertial');
xm1 = [pos*1e-3;vel*1e-3];
%% initialize vars
load('RadarInfo.mat') % sent

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
y_exp(:,1) = gen_observation_fn(xode(:,1),lla,tab.time(1));
%% set up loop
for i = 2:N
    %% sgp4 update
    [pos,vel] = states(sat,tab.time(i),'CoordinateFrame','inertial');
    x_sgp4(:,i) = [pos*1e-3;vel*1e-3];
    
    dt = seconds(tab.time(i)-tab.time(i-1)); 
    Time(i) = seconds(tab.time(i)-tab.time(1));
    Phi_flat = reshape(eye(6),6^2,1);
    z0 = [xode(:,i-1);Phi_flat];
        [~,z] = ode45(@(t,z) dynamics(t,z),[Time(i-1),Time(i)],z0,...
            options);
    xode(:,i) = z(end,1:6);
    
    snsr_idx = find(radar_table.sensor_num==tab.sensor_num(i)); % which radar is this meas from
    lla = [radar_table.Latitude(snsr_idx),radar_table.Longitude(snsr_idx),...
    radar_table.Altitude(snsr_idx)]; % latitude, longitude, and altitude of the radar 
    
    y_exp(:,i) = gen_observation_fn(xode(:,i),lla,tab.time(i));
    y_real(:,i) = [tab.range(i);tab.range_rate(i);tab.az(i);tab.el(i)];
end

%% result 
diff = xode-x_sgp4;
diffmeas = y_real - y_exp;

%% plot
figure('Name','UpdateScriptvODE45','Position',[50,50,700,700])
subplot(3,2,1)
plot(tab.time(1:N),x_sgp4(1,:),'.b','MarkerSize',12),hold on, grid on, grid minor
plot(tab.time(1:N),xode(1,:),'.g')
xlim('tight'),ylim('padded')
ylabel('x km')
legend('sgp4','ode45')
subplot(3,2,3)
plot(tab.time(1:N),x_sgp4(2,:),'.b','MarkerSize',12),hold on, grid on, grid minor
plot(tab.time(1:N),xode(2,:),'.g','MarkerSize',12)
xlim('tight'),ylim('padded')
legend('sgp4','ode45')
ylabel('y km')
subplot(3,2,5)
plot(tab.time(1:N),x_sgp4(3,:),'.b','MarkerSize',12),hold on, grid on, grid minor
plot(tab.time(1:N),xode(3,:),'.g','MarkerSize',12)
xlim('tight'),ylim('padded')
ylabel('z km')
legend('sgp4','ode45')
subplot(3,2,2)
plot(tab.time(1:N),x_sgp4(4,:),'.b','MarkerSize',12),hold on, grid on, grid minor
plot(tab.time(1:N),xode(4,:),'.g','MarkerSize',12)
xlim('tight'),ylim('padded')
ylabel('vx km/s')
legend('sgp4','ode45')
subplot(3,2,4)
plot(tab.time(1:N),x_sgp4(5,:),'.b','MarkerSize',12),hold on, grid on, grid minor
plot(tab.time(1:N),xode(5,:),'.g','MarkerSize',12)
xlim('tight'),ylim('padded')
ylabel('vy km/s')
legend('sgp4','ode45')
subplot(3,2,6)
plot(tab.time(1:N),x_sgp4(6,:),'.b','MarkerSize',12),hold on, grid on, grid minor
plot(tab.time(1:N),xode(6,:),'.g','MarkerSize',12)
xlim('tight'),ylim('padded')
legend('sgp4','ode45')
ylabel('vz km/s')
sgtitle('SGP4 states and ODE45 States')

figure('Name','ODE45tvSGP4diffs','Position',[50,50,700,700])
subplot(3,2,1)
plot(tab.time(1:N),x_sgp4(1,:)-xode(1,:),'.r','MarkerSize',12),hold on, grid on, grid minor
xlim('tight'),ylim('padded')
ylabel('x km')
subplot(3,2,3)
plot(tab.time(1:N),x_sgp4(2,:)-xode(2,:),'.r','MarkerSize',12),hold on, grid on, grid minor
xlim('tight'),ylim('padded')
ylabel('y km')
subplot(3,2,5)
plot(tab.time(1:N),x_sgp4(3,:)-xode(3,:),'.r','MarkerSize',12),hold on, grid on, grid minor
xlim('tight'),ylim('padded')
ylabel('z km')
subplot(3,2,2)
plot(tab.time(1:N),x_sgp4(4,:)-xode(4,:),'.r','MarkerSize',12),hold on, grid on, grid minor
xlim('tight'),ylim('padded')
ylabel('vx km/s')
subplot(3,2,4)
plot(tab.time(1:N),x_sgp4(5,:)-xode(5,:),'.r','MarkerSize',12),hold on, grid on, grid minor
xlim('tight'),ylim('padded')
ylabel('vy km/s')
subplot(3,2,6)
plot(tab.time(1:N),x_sgp4(6,:)-xode(6,:),'.r','MarkerSize',12),hold on, grid on, grid minor
xlim('tight'),ylim('padded')
ylabel('vz km/s')
sgtitle('SGP4 state - ODE45 state')

figure('Name','ExpectedMeasRealMeas','Position',[50,50,700,700])
subplot(4,1,1)
plot(tab.time(1:N),y_real(1,:),'.b','MarkerSize',12),hold on, grid on, grid minor
plot(tab.time(1:N),y_exp(1,:),'.m','MarkerSize',12)
ylabel('Range km')
legend('Real Meas','Expected Meas')
subplot(4,1,2)
plot(tab.time(1:N),y_real(2,:),'.b','MarkerSize',12),hold on, grid on, grid minor
plot(tab.time(1:N),y_exp(2,:),'.m','MarkerSize',12)
xlim('tight'),ylim('padded')
legend('Real Meas','Expected Meas')
ylabel('range-rate km./s')
subplot(4,1,3)
plot(tab.time(1:N),y_real(3,:),'.b','MarkerSize',12),hold on, grid on, grid minor
plot(tab.time(1:N),y_exp(3,:),'.m','MarkerSize',12)
xlim('tight'),ylim('padded')
ylabel('el, deg')
legend('Real Meas','Expected Meas')
subplot(4,1,4)
plot(tab.time(1:N),y_real(4,:),'.b','MarkerSize',12),hold on, grid on, grid minor
plot(tab.time(1:N),y_exp(4,:),'.m','MarkerSize',12)
xlim('tight'),ylim('padded')
ylabel('az, deg')
legend('Real Meas','Expected Meas')
sgtitle('Real Measurements and Expected Measurements')

figure('Name','ExpectedMeasRealMeasDIFF','Position',[50,50,700,700])
subplot(4,1,1)
plot(tab.time(1:N),diffmeas(1,:),'.r','MarkerSize',15),hold on, grid on, grid minor
xlim('tight'),ylim('padded')
ylabel('Range km')
subplot(4,1,2)
plot(tab.time(1:N),diffmeas(2,:),'.r','MarkerSize',15),hold on, grid on, grid minor
xlim('tight'),ylim('padded')
ylabel('range-rate km./s')
subplot(4,1,3)
plot(tab.time(1:N),diffmeas(3,:),'.r','MarkerSize',15),hold on, grid on, grid minor
xlim('tight'),ylim('padded')
ylabel('el, deg')
subplot(4,1,4)
plot(tab.time(1:N),diffmeas(4,:),'.r','MarkerSize',15),hold on, grid on, grid minor
xlim('tight'),ylim('padded')
ylabel('az, deg')
sgtitle('Real Measurements minus Expected Measurements')

