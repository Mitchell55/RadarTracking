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

%% Select Satellite
tab = tab(tab.fin_tag == 47492,:); % change as needed

%% get the  day frac to set up SGP4
% SGP-4 propagator is used to propagate TLEs. TLEs need to have fractional
% days
Tfrac = tab.day(1) + tab.hr(1)/24 + tab.min(1)/1440 + tab.s(1)/86400;
year = tab.year(1);

t1 = datetime(tab.year(1) + 2000,1,tab.day(1),tab.hr(1),tab.min(1),tab.s(1));
t2 = datetime(tab.year(2) + 2000,1,tab.day(2),tab.hr(2),tab.min(2),tab.s(2));

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
    seconds(t2 - t1), 0, 398600.435507);

% find closest velocity to 7 km/sec -- should just be v2
[val, idx] = min(abs([norm(V1), norm(V2)] - 7));
vel = V2;

% dt is roughly 10 seconds - radar meases for rougly 170 seconds
state = [eci_pos_2; vel'];
Phi_flat = reshape(eye(6),6^2,1);
state = [state;Phi_flat];
t_sim = [datetime(datevec(t2))];
% range, range rate, az (deg), el (deg)
sig_range = 5e-3;
sig_range_rate = 1e-3;
sig_ang = 1e-3;
% Process noise matrix
dt = 10;
sig_q = 1e-6;
Q = ((sig_q/10)^2)*...
    [((dt^3)/3), 0, 0, ((dt^2)/2), 0, 0;
    0, ((dt^3)/3), 0, 0, ((dt^2)/2), 0;
    0, 0, ((dt^3)/3), 0, 0, ((dt^2)/2);
    ((dt^2)/2), 0, 0, dt, 0, 0;
    0, ((dt^2)/2), 0, 0, dt, 0;
    0, 0, ((dt^2)/2), 0, 0, dt];
% add large offset
radar_observation = gnc.gen_observation_fn(state(1:6,end) + [0, 0, 0, 0, 0, 0]', lla2, t2);
for t = 10:10:200
    [prop_times,final_state] = ode45(@(t,final_state) gnc.state_dyn(t,final_state),...
        [0, 10],state(:,end),...
        options);
    t_sim = [t_sim; t_sim(1) + seconds(t)];
    state = [state, final_state(end,:)'];
    state(7:end,end) = Phi_flat;
    process_noise = mvnrnd([0,0,0,0,0,0],Q)';
    state(1:6,end) = state(1:6,end) + process_noise;
    radar_observation = [radar_observation, ...
        gnc.gen_observation_fn(state(1:6,end),lla2, t_sim(end) )];
    % add measurement noise
    radar_observation(1,end) = radar_observation(1,end) + 0*sig_range*randn;
    radar_observation(2,end) = radar_observation(2,end) + 0*sig_range_rate*randn;
    radar_observation(3,end) = radar_observation(3,end) + 0*sig_ang*randn;
    radar_observation(4,end) = radar_observation(4,end) + 0*sig_ang*randn;
end

configs.R = diag( [ (10e-3)^2, (10e-3)^2, (10e-3)^2] );
filter = batchFilter.BatchFilter(configs);

% [initial_state, filter_residuals, eciMeas] = filter.process_meas(radar_observation,...
%     t_sim, repmat(lla2,[size(radar_observation, 2), 1] ), vel' )

% init the filter at the wrong initial state
filtered_states = state;
filtered_states(1:6,1) = filtered_states(1:6,1) + [5, 0, 0, 0, 0, 0]';

initial_state = filter.getInitialState( filtered_states(1:3,:), t_sim, vel' );
filter_residuals = filter.getFilterResiduals( initial_state, state(1:3,:), t_sim, vel' );


%% plot the filter residuals
figure;
subplot(1,3,1)
plot(t_sim,filter_residuals(1,:), '.')
title('x')
ylabel('m')
xlabel('s')

subplot(1,3,2)
plot(t_sim,filter_residuals(2,:), '.')
title('y')
ylabel('m')
xlabel('s')

subplot(1,3,3)
plot(t_sim,filter_residuals(3,:), '.')
title('z')
ylabel('m')
xlabel('s')


%% reset - runnin an additional 8 times seems to make it work
for i = 1:1
filtered_states(1:6,1) = initial_state;

initial_state = filter.getInitialState( filtered_states(1:3,:), t_sim, vel' );
filter_residuals = filter.getFilterResiduals( initial_state, state(1:3,:), t_sim, vel' );
end

%% Experiment - Interesting

tmp_state = state(1:6,3);
radar_obs = gnc.gen_observation_fn(tmp_state,lla2, t_sim(1) );
convert_state = gnc.getstate(radar_obs,lla2,t_sim(1))

error = tmp_state(1:3) - convert_state
norm(error)*1e3 % off by 47 meters?


