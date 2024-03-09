function [r_eci, r_ecef] = getstate(Z,lla,time)
az = Z(3); % deg
el = Z(4);
% range_rate = Y(2);
rho = Z(1);

rho_sez(1) = -rho*cosd(el)*cosd(az);
rho_sez(2) = rho*cosd(el)*sind(az);
rho_sez(3) = rho*sind(el);

% go from SEZ to ECI
C_sez2ecef = [sind(lla(1))*cosd(lla(2)), -sind(lla(2)), cosd(lla(1))*cosd(lla(2));
    sind(lla(1))*sind(lla(2)), cosd(lla(2)), cosd(lla(1))*sind(lla(2));
    -cosd(lla(1)), 0, sind(lla(1))];

rho_ecef = C_sez2ecef*rho_sez';
% now get site ecef
% all units km
% R_earth = 6378.137; %km
% e_earth = 0.081819;
% h_ell = (lla(3) + egm96geoid(lla(1),lla(2)))*1e-3;
% C_earth = R_earth/sqrt(1-e_earth^2*(sind(lla(1)))^2);
% S_earth = C_earth*(1-e_earth^2);
% 
% r_delt = (C_earth+h_ell)*cosd(lla(1));
% r_K = (S_earth+h_ell)*sind(lla(1));
% site_ecef = [r_delt*cosd(lla(2));r_delt*sind(lla(2));r_K];
site_ecef =  lla2ecef(lla)'/1e3;

% get pos in ecef;
X_ecef = rho_ecef+site_ecef;
r_ecef = X_ecef;

% go to ECI
[r_eci,~,~] = ecef2eci(datevec(time),X_ecef);

end
