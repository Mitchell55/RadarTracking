load("data\obs_split.mat")

idx = 1;

Z = obs_table{idx,:};
sensor = obs_table{idx,2};
sensor_idx = (radar_table{:,2}==sensor);
lat = radar_table{sensor_idx, 3};
lon = radar_table{sensor_idx, 4};
alt = radar_table{sensor_idx, 5};
lla = [lat, lon, alt];

year = obs_table{idx,3};
day = obs_table{idx, 4};
hr = obs_table{idx, 5};
minu = obs_table{idx, 6};
s = obs_table{idx, 7};
t = datenum(year + 2000,0,day,hr,minu,s);

[eci, ecef] = gnc.getstate(Z, lla,t);

options = odeset('RelTol',1e-13,'AbsTol',1e-15);
[~,z] = ode45(@(t,z) gnc.state_dyn(t,z),[t,t + 100],eci,...
            options);
 



