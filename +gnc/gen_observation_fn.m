function Z= gen_observation_fn(X,lla,time)

[R,V,~] = eci2ecef(datevec(time),X(1:3).*1e3,X(4:6).*1e3); % you were here
R = R*1e-3; V = V*1e-3;

Rs = lla2ecef(lla);
Rs = Rs.*1e-3;

range_vec = R-Rs';
range = norm(range_vec);
rngrt_vec = V;

latitude = lla(1)*pi/180;
longitude = lla(2)*pi/180;

C_lat = [cos(pi/2-latitude),0,-sin(pi/2-latitude);
    0,1,0;
    sin(pi/2-latitude),0,cos(pi/2-latitude)];
C_lon = [cos(longitude),sin(longitude),0;
    -sin(longitude),cos(longitude),0;
    0,0,1];
C_ecef2sez = C_lat*C_lon;

rho_sez = C_ecef2sez*range_vec;
rho_sez_dot = C_ecef2sez * rngrt_vec;

range_rate = dot(rho_sez,rho_sez_dot,1)./range;
el = asin(rho_sez(3)/range);
az = atan2(rho_sez(2) , -rho_sez(1) ); % y,x_ecef - cos,sin
az(az<0) = 2*pi+az;
Z = [range;range_rate;az*180/pi;el*180/pi];


end