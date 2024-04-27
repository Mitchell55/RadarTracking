%% Load data
load("data\obs_split.mat")

%% calculate ECI location
meas_idx = 1;
sensor_num = obs_table{meas_idx, 2};
sensor_idx = (radar_table{:,2}==sensor_num);
lat = radar_table{sensor_idx, 3};
lon = radar_table{sensor_idx, 4};
alt = radar_table{sensor_idx, 5};

Z = obs_table{meas_idx,:};
lla = [lat, lon, alt];
year = obs_table{meas_idx,3};
day = obs_table{meas_idx, 4};
hr = obs_table{meas_idx, 5};
minu = obs_table{meas_idx, 6};
s = obs_table{meas_idx, 7};
t = datenum(year + 2000,0,day,hr,minu,s);

[r v] = gnc.getstate(Z,lla,datetime([year(1) + 2000, 1, day(1), hr(1), minu(1), s(1)]))




%% trying to find overlapping regions
% count = 0;
% times = cell(1);
% sensor_list = cell(1);
% for target = unique(obs_table{:,1})'
%     for sensor = unique(obs_table{:,2})'
%         idx = find(obs_table{:,1}==target & obs_table{:,2}==sensor);
%         if ~isempty(idx)
%             count = count+1;
%             year = obs_table{idx,3};
%             day = obs_table{idx, 4};
%             hr = obs_table{idx, 5};
%             minu = obs_table{idx, 6};
%             s = obs_table{idx, 7};
%             t = datenum(year + 2000,0,day,hr,minu,s);
%             sensor_num = sensor*ones(size(t));
%             times{count} =  t;
%             sensor_list{count} =  sensor_num';
%             
%         end
%         
%     end
%     % now plot the data
%     figure;
%     hold on
%     for i = 1:numel(sensor_list)
%         plot(times{i},sensor_list{i},'x')
%     end
%     title(['target ' num2str(target)] )
% end

%% plot the eci data for target 1

target = obs_table{1,1};
count = 0;
times = cell(1);
sensor_list = cell(1);
eci_x = cell(1);
eci_y = cell(1);
eci_z = cell(1);
lat_data = cell(1);
lon_data = cell(1);
alt_data = cell(1);

for sensor = unique(obs_table{:,2})'
    idx = (obs_table{:,1}==target & obs_table{:,2}==sensor);
    if sum(idx~=0)
        year = obs_table{idx,3};
        day = obs_table{idx, 4};
        hr = obs_table{idx, 5};
        minu = obs_table{idx, 6};
        s = obs_table{idx, 7};
        t = datenum(year + 2000,0,day,hr,minu,s);
        
        sensor_idx = (radar_table{:,2}==sensor);
        lat = radar_table{sensor_idx, 3};
        lon = radar_table{sensor_idx, 4};
        alt = radar_table{sensor_idx, 5};
        
        Z = obs_table{idx,:};
        lla = [lat, lon, alt];
        if ~isnan(lla)
            count = count + 1;
            eci = [];
            ecef = [];
            for i = 1:size(t,1)
                [eci(:,i), ecef(:,i)] = gnc.getstate(Z(i,:),lla,t(i));
            end
            lla_meas = ecef2lla(ecef'*1e3);
            times{count} =  t;
            sensor_list{count} =  sensor;
            eci_x{count} = eci(1,:);
            eci_y{count} = eci(2,:);
            eci_z{count} = eci(3,:);
            lat_data{count} = lla_meas(:,1);
            lon_data{count} = lla_meas(:,2);
            alt_data{count} = lla_meas(:,3);
            count
        end
    end
end



minTime = inf;
for k = 1:numel(times)
    min_t = min(times{k});
    if min_t<minTime
        minTime = min_t;
    end
end


%% total data

figure;
leg = {};
hold on
for k = 1:numel(times)
    leg{k} = num2str(sensor_list{k})
    plot((times{k} - minTime)*86400,lat_data{k},'x','linewidth',10)
end
legend(leg)
hold off
grid on
grid minor
title(['satellite ' num2str(target)]);
ylabel('Lattitude')
xlabel('Days')



%% ground track
figure;
leg = {};
hold on
for k = 1:numel(times)
    leg{k} = num2str(sensor_list{k})
    t = (times{k} - minTime)*86400;
    lat = lat_data{k};
    lon = lon_data{k};
    a = alt_data{k};
%     idx = (t>1.1025e4 +90*60 | t< 1.1025e4);
%     t(idx) = [];
%     lat(idx) = [];
%     lon(idx) = [];
      plot(lon_data{k}, lat_data{k})
%     plot(t,lat_data{k},'x','linewidth',6)
%     plot(a,'x','linewidth',5)
end
legend(leg)
hold off
grid on
grid minor
title(['satellite ' num2str(target)]);
ylabel('Lattitude')
xlabel('Longitude')

%% Try to propagate
x = eci_x{1};
y = eci_y{1};
z = eci_z{1};

t_start = times{1}*86400;

pos = [x(end), y(end), z(end)];
vel = [x(end) - x(end-1), y(end) - y(end-1), ...
    z(end) - z(end-1)]./(t_start(end) - t_start(end - 1));

% state = [pos'; vel'];
state = [1024.17896349139;
-3219.32174203467;
6017.86991764859;
-0.0320477794777778;
-6.70180835785332;
-3.58085289433470];

Phi_flat = reshape(eye(6),6^2,1);
state = [state;Phi_flat];
tf = times{3}*86400;
tf = tf(end);

options = odeset('RelTol',1e-13,'AbsTol',1e-15);

[prop_times,final_state] = ode45(@(t,final_state) gnc.state_dyn(t,final_state),...
    [0,90*60],state,...
            options);

final_eci = [final_state(:,1), final_state(:,2),... 
    final_state(:,3)];

prop_ecf = [];
for i = 1:100:size(final_eci,1)
    prop_ecf = [prop_ecf; eci2ecef(datevec(t_start(end)/86400 + prop_times(i)/86400),...
        final_eci(i,:).*1e3)'];
end

prop_lla = ecef2lla(prop_ecf);

figure; plot(prop_lla(:,2), prop_lla(:,1),'.')
figure; plot(prop_times, final_state(:,1),'.')
%%
figure;
leg = {};
hold on
count = 0;
for k = 1:numel(times)
    t = (times{k})*86400;
    idx = ( ( t>t_start(end))&(t < t_start(end) + 90*60) );
    if sum(idx>0)
        count = count + 1;
        leg{count} = num2str(sensor_list{k})
        plot(lon_data{k}, lat_data{k})
    end
end
plot(prop_lla(:,2), prop_lla(:,1),'--')
legend(leg)
hold off
grid on
grid minor
title(['satellite ' num2str(target)]);
ylabel('Lattitude')
xlabel('Longitude')

