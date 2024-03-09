%% Clearing
clc; clear, close all

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
% tab = tabp{400};

%% Select Satellite
tab = tab(tab.fin_tag == 47492,:); % change as needed 47422 - starlink
% tab = tab(tab.fin_tag == 47462, :);

%% get the  day frac to set up SGP4
% SGP-4 propagator is used to propagate TLEs. TLEs need to have fractional
% days
Tfrac = tab.day(1) + tab.hr(1)/24 + tab.min(1)/1440 + tab.s(1)/86400;
year = tab.year(1);

t1 = datenum(tab.year(1) + 2000,1,tab.day(1),tab.hr(1),tab.min(1),tab.s(1));
t2 = datenum(tab.year(2) + 2000,1,tab.day(2),tab.hr(2),tab.min(2),tab.s(2));

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


[eci_pos_1, ecef_pos_1] = gnc.getstate(z1,lla1,t1);
[eci_pos_2, ecef_pos_2] = gnc.getstate(z2,lla2,t2);

% now solve for initial state using lambert
[V1, V2, extremal_distances, exitflag] = gnc.lambert(eci_pos_1', eci_pos_2',...
    t2 - t1, 0, 398600.435507);

% find closest velocity to 7 km/sec -- should just be v2
[val, idx] = min(abs([norm(V1), norm(V2)] - 7));
vel = V2;

% states, times, monte carlo run
m = 1;

% dt is roughly 10 seconds - radar meases for rougly 170 seconds
start_state = [eci_pos_2; vel'];
Phi_flat = reshape(eye(6),6^2,1);
start_state = [start_state;Phi_flat];
state = nan(42, 20, m);
for i = 1:m
    state(:,1,i) = start_state;
end

t_sim = [datetime(datevec(t2))];
% range, range rate, az (deg), el (deg)
sig_range = 5e-3;
sig_range_rate = 1e-3;
sig_ang = 1e-3;
% Process noise matrix
dt = 10; %10
sig_q = 1e-6;
Q = ((sig_q/10)^2)*...
    [((dt^3)/3), 0, 0, ((dt^2)/2), 0, 0;
    0, ((dt^3)/3), 0, 0, ((dt^2)/2), 0;
    0, 0, ((dt^3)/3), 0, 0, ((dt^2)/2);
    ((dt^2)/2), 0, 0, dt, 0, 0;
    0, ((dt^2)/2), 0, 0, dt, 0;
    0, 0, ((dt^2)/2), 0, 0, dt];

radar_observation = nan(4,20, m);
for monte = 1:m
    count = 1;
    radar_observation(:,1,monte) = gnc.gen_observation_fn(state(1:6, 1,monte),lla2, t2);
    for t = dt:dt:200
        count = count + 1;
        [prop_times,final_state] = ode45(@(t,final_state) gnc.state_dyn(t,final_state),...
            [0, dt],state(:,count-1, monte),...
            options);
        t_sim = [t_sim; t_sim(1) + seconds(t)];
        state(:,count, monte) = final_state(end,:)';
        state(7:end,count, monte) = Phi_flat;
        process_noise = mvnrnd([0,0,0,0,0,0],Q)';
        state(1:6,count,end) = state(1:6,count,end) + 0*process_noise;
        radar_observation(:,count, monte) = gnc.gen_observation_fn(state(1:6,count, monte), lla2, t_sim(end));
        % add measurement noise
        radar_observation(1,end) = radar_observation(1,end) + 0*sig_range*randn;
        radar_observation(2,end) = radar_observation(2,end) + 0*sig_range_rate*randn;
        radar_observation(3,end) = radar_observation(3,end) + 0*sig_ang*randn;
        radar_observation(4,end) = radar_observation(4,end) + 0*sig_ang*randn;
    end
end

% 
% final_eci = [state(1,:)', state(2,:)',... 
%     state(3,:)'];
% 
% prop_ecf = [];
% for i = 1:size(final_eci,1)
%     prop_ecf = [prop_ecf; eci2ecef(datevec(t_sim(i)),...
%         final_eci(i,:).*1e3)'];
% 
% end


%% Run the EKF
ekf_state = nan(6,20,m);
state_error = ekf_state;
filter_resid = nan(4,20,m);
P = nan(6,6, size(radar_observation,2), m );
sigmas = nan(6, size(radar_observation, 2 ), monte );
for monte = 1:m
    clear configs
    clear filter
    configs = configs.sim_config(state(1:6,1,monte) + [0; 0; 5; 0; 0; 0], sig_q);
    % configs = configs.sim_config(state(1:6,1,monte) + [0.001; -0.008; -0.009; 0; 0; 0], sig_q);
    filter = ekf.EKF(configs, t_sim(1));
    
    for i = 2:size(radar_observation,2)
        resid = filter.process_meas( radar_observation(:,i,monte), t_sim(i), ...
            lla2 );
        filter_resid(:,i,monte) = resid;
        ekf_state(:,i,monte) = filter.state;
        state_error(:,i, monte) = filter.state - state(1:6,i,monte);
        P(:,:, i, monte) = filter.P;
        sigmas(:,i, monte) = diag(P(:,:,i,monte)).^(1/2);
    end
end




%% Plot the error in the filtering 

average_error = sum(state_error,3)/m;
average_sigmas = sum(sigmas, 3)/m;

figure;
subplot(2,3,1)
plot(average_error(1,:),'.b', 'MarkerSize',9)
hold on
plot(average_sigmas(1,:),'r--')
plot(-average_sigmas(1,:),'r--')
hold off
grid on 
grid minor
ylabel('Km')
title('x error')
legend('State error', '1 sigma bounds')

subplot(2,3,2)
plot(average_error(2,:),'.b', 'MarkerSize',9)
hold on
plot(average_sigmas(2,:),'r--')
plot(-average_sigmas(2,:),'r--')
hold off
grid on 
grid minor
title('y error')
ylabel('Km')

subplot(2,3,3)
plot(average_error(3,:),'.b', 'MarkerSize',9)
hold on
plot(average_sigmas(3,:),'r--')
plot(-average_sigmas(3,:),'r--')
hold off
grid on 
grid minor
title('z error')
ylabel('Km')

subplot(2,3,4)
plot(average_error(4,:),'.b', 'MarkerSize',9)
hold on
plot(average_sigmas(4,:),'r--')
plot(-average_sigmas(4,:),'r--')
hold off
grid on 
grid minor
title('x vel error')
ylabel('Km/s')

subplot(2,3,5)
plot(average_error(5,:),'.b', 'MarkerSize',9)
hold on
plot(average_sigmas(5,:),'r--')
plot(-average_sigmas(5,:),'r--')
hold off
grid on 
grid minor
title('y vel error')
ylabel('Km/s')

subplot(2,3,6)
plot(average_error(6,:),'.b', 'MarkerSize',9)
hold on
plot(average_sigmas(6,:),'r--')
plot(-average_sigmas(6,:),'r--')
hold off
grid on 
grid minor
title('z vel error')
ylabel('Km/s')

%% Plot the average filter residuals 

average_resid = sum(filter_resid,3)/m;
average_sigmas = sum(sigmas, 3)/m;

figure;
subplot(2,2,1)
plot(average_resid(1,:),'.b', 'MarkerSize',9)
grid on 
grid minor
ylabel('Km')
title('range residuals')

subplot(2,2,2)
plot(average_resid(2,:),'.b', 'MarkerSize',9)
grid on 
grid minor
title('range rate residuals')
ylabel('Km/s')

subplot(2,2,3)
plot(average_resid(3,:),'.b', 'MarkerSize',9)
grid on 
grid minor
title('az residuals')
ylabel('deg')

subplot(2,2,4)
plot(average_resid(4,:),'.b', 'MarkerSize',9)
grid on 
grid minor
title('el residuals')
ylabel('deg')
