% filter actual radar data
clc
clear all
close all


%% Figure settings
set(0,'DefaultFigureColor',[1 1 1]); set(0,'DefaultLineLineWidth',0.9);
set(0,'DefaultAxesFontSize',12); set(0,'DefaultTextInterpreter','latex');
set(0,'DefaultAxesFontWeight','Normal');
set(0,'DefaultAxesTitleFontWeight','Bold');
%% load in info
load('data\obs_split')

%% ODE settings
options = odeset('RelTol',1e-13,'AbsTol',1e-15);

%% Select Pass
tab = tabp{300}; % looking at pass 300 - Fylingdales SW

%% Select Satellite
tab = tab(tab.fin_tag == 47492,:); % change as needed

%% get all the times and make sure they are sorted
% SGP-4 propagator is used to propagate TLEs. TLEs need to have fractional
% days
Tfrac = tab.day(:) + tab.hr(:)/24 + tab.min(:)/1440 + tab.s(:)/86400;
year = tab.year(:);

% looks like times are already sorted
times = datetime(tab.year(:) + 2000,1,tab.day(:),tab.hr(:),tab.min(:),tab.s(:));

%% solve for the initial state
sensor_num_1 = tab.sensor_num(1);
sensor_num_2 = tab.sensor_num(2);

sensor_idx_1 = radar_table.sensor_num==sensor_num_1;
sensor_idx_2 = radar_table.sensor_num==sensor_num_2;

lla1 = [radar_table.Latitude(sensor_idx_1), radar_table.Longitude(sensor_idx_1),...
    radar_table.Altitude(sensor_idx_1)];
lla2 = [radar_table.Latitude(sensor_idx_2), radar_table.Longitude(sensor_idx_2),...
    radar_table.Altitude(sensor_idx_2)];


% tle_file = 'TLEs\47492_12.tle';
% [sc, sat, sat_ID] = gnc.generate_sgp4_scenario(tab.time(1),tab.time(end),...
%     1,tle_file,'sgp4');

% from what I can tell, Z is range, range rate, az, and el
z1 = zeros(1,4);
z2 = z1;
z1(1) = tab.range(1);
z1(2) = tab.range_rate(1);
z1(3) = tab.az(1);
z1(4) = tab.el(1);

z2(1) = tab.range(2);
z2(2) = tab.range_rate(2);
z2(3) = tab.az(2);
z2(4) = tab.el(2);


[eci_pos_1, ecef_pos_1] = gnc.getstate(z1,lla1,times(1));
[eci_pos_2, ecef_pos_2] = gnc.getstate(z2,lla2,times(2) );

% now solve for initial state using lambert
[V1, V2, extremal_distances, exitflag] = gnc.lambert(eci_pos_1', eci_pos_2',...
    datenum(times(2) - times(1) ), 0, 398600.435507);

% find closest velocity to 7 km/sec -- should just be v2
[val, idx] = min(abs([norm(V1), norm(V2)] - 7));
vel = V2;

%% init the filter
start_state = [eci_pos_2; vel'];
% Phi_flat = reshape(eye(6),6^2,1);
start_state = [start_state];
state = nan(6, size(tab,1) - 1);
P = nan(6,6,size(tab,1) - 1);
sigmas = nan(6, size(tab,1) - 1 );
state(:,1) = start_state;

clear configs;
clear filter;

sig_q = 1e-6;
configs = configs.real_data_config(state(1:6,1), sig_q);
P(:,:,1) = configs.P;
sigmas(:,1) = diag(P(:,:,1)).^(1/2);
filter = ekf.EKF(configs, times(2));
radar_observation = nan(4,1);
filter_resid = nan(4,size(tab,1) - 2);
for i = 2:size(tab,1)-1
    % get radar observation from data
    el = tab{i+1, 8};
    az = tab{i+1, 9};
    range = tab{i+1, 10};
    range_rate = tab{i+1, 11};
    radar_observation(1) = range; radar_observation(2) = range_rate;
    radar_observation(3) = az; radar_observation(4) = el;
    % get time from data
    new_time = datetime(tab.year(i+1) + 2000,1,tab.day(i+1),tab.hr(i+1),tab.min(i+1),tab.s(i+1));
    % get the location of the radar
    sensor_num = tab.sensor_num(i + 1);
    sensor_idx = radar_table.sensor_num==sensor_num;
    lla = [radar_table.Latitude(sensor_idx), radar_table.Longitude(sensor_idx),...
        radar_table.Altitude(sensor_idx)];    

    % filter data
    resid = filter.process_meas( radar_observation, new_time, ...
            lla);
        filter_resid(:,i-1) = resid;
        state(:,i) = filter.state;
        P(:,:, i) = filter.P;
        sigmas(:,i) = diag(P(:,:,i)).^(1/2);
end

%% Plot the average filter residuals 

average_resid = filter_resid;
average_sigmas = sigmas;

figure;
subplot(2,2,1)
plot(times(3:end,:),average_resid(1,:),'.b', 'MarkerSize',9)
grid on 
grid minor
ylabel('Km')
title('range residuals')

subplot(2,2,2)
plot(times(3:end,:),average_resid(2,:),'.b', 'MarkerSize',9)
grid on 
grid minor
title('range rate residuals')
ylabel('Km/s')

subplot(2,2,3)
plot(times(3:end,:),average_resid(3,:),'.b', 'MarkerSize',9)
grid on 
grid minor
title('az residuals')
ylabel('deg')

subplot(2,2,4)
plot(times(3:end,:),average_resid(4,:),'.b', 'MarkerSize',9)
grid on 
grid minor
title('el residuals')
ylabel('deg')


%% Plot the states

figure;
subplot(2,3,1)
plot(times(2:end,:),state(1,:),'.b', 'MarkerSize',9)
hold on
plot(times(2:end,:),state(1,:) + average_sigmas(1,:),'r--')
plot(times(2:end,:), state(1,:) -average_sigmas(1,:),'r--')
hold off
grid on 
grid minor
ylabel('Km')
title('x')
legend('State error', '1 sigma bounds')

subplot(2,3,2)
plot(times(2:end,:),state(2,:),'.b', 'MarkerSize',9)
hold on
plot(times(2:end,:),state(2,:) + average_sigmas(2,:),'r--')
plot(times(2:end,:),state(2,:) -average_sigmas(2,:),'r--')
hold off
grid on 
grid minor
title('y')
ylabel('Km')

subplot(2,3,3)
plot(times(2:end,:),state(3,:),'.b', 'MarkerSize',9)
hold on
plot(times(2:end,:),state(3,:) + average_sigmas(3,:),'r--')
plot(times(2:end,:),state(3,:) -average_sigmas(3,:),'r--')
hold off
grid on 
grid minor
title('z')
ylabel('Km')

subplot(2,3,4)
plot(times(2:end,:),state(4,:),'.b', 'MarkerSize',9)
hold on
plot(times(2:end,:),state(4,:) + average_sigmas(4,:),'r--')
plot(times(2:end,:),state(4,:) -average_sigmas(4,:),'r--')
hold off
grid on 
grid minor
title('x vel')
ylabel('Km/s')

subplot(2,3,5)
plot(times(2:end,:),state(5,:),'.b', 'MarkerSize',9)
hold on
plot(times(2:end,:),state(5,:) + average_sigmas(5,:),'r--')
plot(times(2:end,:),state(5,:) -average_sigmas(5,:),'r--')
hold off
grid on 
grid minor
title('y vel')
ylabel('Km/s')

subplot(2,3,6)
plot(times(2:end,:),state(6,:),'.b', 'MarkerSize',9)
hold on
plot(times(2:end,:),state(6,:) + average_sigmas(6,:),'r--')
plot(times(2:end,:),state(6,:) -average_sigmas(6,:),'r--')
hold off
grid on 
grid minor
title('z vel')
ylabel('Km/s')

%% measurement sub plot
% 
% figure;
% subplot(2,2,1)
% plot(times(2:end,:),,'.b', 'MarkerSize',9)
% hold on
% plot(times(2:end,:),state(1,:) + average_sigmas(1,:),'r--')
% plot(times(2:end,:), state(1,:) -average_sigmas(1,:),'r--')
% hold off
% grid on 
% grid minor
% ylabel('Km')
% title('x')
% legend('State error', '1 sigma bounds')
% 
% subplot(2,3,2)
% plot(times(2:end,:),state(2,:),'.b', 'MarkerSize',9)
% hold on
% plot(times(2:end,:),state(2,:) + average_sigmas(2,:),'r--')
% plot(times(2:end,:),state(2,:) -average_sigmas(2,:),'r--')
% hold off
% grid on 
% grid minor
% title('y')
% ylabel('Km')
% 
% subplot(2,3,3)
% plot(times(2:end,:),state(3,:),'.b', 'MarkerSize',9)
% hold on
% plot(times(2:end,:),state(3,:) + average_sigmas(3,:),'r--')
% plot(times(2:end,:),state(3,:) -average_sigmas(3,:),'r--')
% hold off
% grid on 
% grid minor
% title('z')
% ylabel('Km')
% 
% subplot(2,3,4)
% plot(times(2:end,:),state(4,:),'.b', 'MarkerSize',9)
% hold on
% plot(times(2:end,:),state(4,:) + average_sigmas(4,:),'r--')
% plot(times(2:end,:),state(4,:) -average_sigmas(4,:),'r--')
% hold off
% grid on 
% grid minor
% title('x vel')
% ylabel('Km/s')
% 
% subplot(2,3,5)
% plot(times(2:end,:),state(5,:),'.b', 'MarkerSize',9)
% hold on
% plot(times(2:end,:),state(5,:) + average_sigmas(5,:),'r--')
% plot(times(2:end,:),state(5,:) -average_sigmas(5,:),'r--')
% hold off
% grid on 
% grid minor
% title('y vel')
% ylabel('Km/s')
% 
% subplot(2,3,6)
% plot(times(2:end,:),state(6,:),'.b', 'MarkerSize',9)
% hold on
% plot(times(2:end,:),state(6,:) + average_sigmas(6,:),'r--')
% plot(times(2:end,:),state(6,:) -average_sigmas(6,:),'r--')
% hold off
% grid on 
% grid minor
% title('z vel')
% ylabel('Km/s')
% 
