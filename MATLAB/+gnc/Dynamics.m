%% Get dynamics of system
%% clearing
clear
close all

%% symbolic toolbox
syms J2 x y z mu dx dy dz rE real
state = [x,y,z,dx,dy,dz];

% converting to cartesian
sinphi = z/sqrt(x^2+y^2+z^2);
% calculate the coefficients
P0 = 1/sinphi;
P1 = 1;
P2 = (3*sinphi*P1*sinphi-1*P0*sinphi)/(2*sinphi);
r = sqrt(x^2+y^2+z^2); % define orbital radius
U(x,y,z) = (mu/r)*(1-(rE/r)^2*P2*sinphi*J2); % potential energy
acc = gradient(U(x,y,z),[x,y,z]); % spacecraft accelerations 

ddx = acc(1);
ddy = acc(2);
ddz = acc(3);
dx = dx;
dy = dy;
dz = dz;
f = [dx,dy,dz,ddx,ddy,ddz];
A = jacobian(f,[x,y,z,dx,dy,dz]);

%% now the acc_drag 
% acc_dragecef = -(1/2)*rho*(Cd*Area/m).*norm(Vr)*Vr; % this is the acceleration in the ecef position.
%         % Need transportation theorem to go back and forth 
% % acc_drag = C_eci2ecef'*acc_dragecef + cross([0,0,wdot]',Vr); % guess on the transportation theorem
% acc_drag = acc_dragecef;
% acc = acc_nodrag + acc_drag;

% for part 1 we are just interested in the 6x1 sc state
% ddx = acc(1);
% ddy = acc(2);
% ddz = acc(3);
% dx = dx;
% dy = dy;
% dz = dz;
% f = [dx,dy,dz,ddx,ddy,ddz];
% A = jacobian(f,[x,y,z,dx,dy,dz]);

%% Part 2: h matrix
syms xs ys zs dxs dys dzs lat_s long_s om_E real
syms Uxx Uxy Uxz Uyx Uyy Uyz Uzx Uzy Uzz real
syms dUxx dUxy dUxz dUyx dUyy dUyz dUzx dUzy dUzz real
R_eci = [x;y;z]; V_eci = [dx;dy;dz];
% U = [Uxx,Uxy,Uxz;
%     Uyx,Uyy,Uyz;
%     Uzx,Uzy,Uzz];
% dU = [dUxx,dUxy,dUxz;
%     dUyx,dUyy,dUyz;
%     dUzx,dUzy,dUzz];
% 
% R = U*R_eci;
% V = U*V_eci + dU*R_eci;
R = R_eci; V = V_eci;
Rs = [xs;ys;zs]; % ecef


range_vec = R-Rs;
rngrt_vec = V;
range = norm(range_vec);

latitude = lat_s; % rad
longitude = long_s; % rad

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
y_ex = [norm(rho_sez);dot(rho_sez,rho_sez_dot)./range;el;az];


% h =[range;range_rate;az*180/pi;el*180/pi];
h =[az*180/pi;el*180/pi]
H = jacobian(h,[x,y,z,dx,dy,dz])

%%

% syms x y z dx dy dz xu yu zu dxu dyu dzu real
% 
% r_ijk_sat = [x;y;z]; v_ijk_sat = [dx;dy;dz];
% r_ijk_user = [xu;yu;zu]; v_ijk_user = [dxu;dyu;dzu];
% 
% rho_ijk = r_ijk_sat - r_ijk_user; rho = norm(rho_ijk);
% rho_dot_ijk = v_ijk_sat-v_ijk_user; rho_dot = dot(rho_ijk,rho_dot_ijk)/rho;
% 
% dec_topo = asin(rho_ijk(3)/rho);
% ra_topo = atan2(rho_ijk(2),rho_ijk(1));
% 
% dec_dot_topo = (rho_dot_ijk(3)-rho_dot*sin(dec_topo))/sqrt(rho_ijk(1)^2 + rho_ijk(2)^2);
% ra_dot_topo = (rho_dot_ijk(1)*rho_ijk(2)-rho_dot_ijk(2)*rho_ijk(2))/(-rho_ijk(2)^2 - rho_ijk(1)^2);
% 
% 
% ymeas = [dec_topo;ra_topo;dec_dot_topo;ra_dot_topo];
% 
% H = jacobian(ymeas,[xu,yu,zu,dxu,dyu,dzu])
% 
% 
% %% messing around
% % H3
% pretty((-(1/((x-xu)^2 + (y-yu)^2 + (z-zu)^2)^(1/2) - ((z - zu)*(z - zu))/((x-xu)^2 + (y-yu)^2 + (z-zu)^2)^(3/2))/(1 - (z - zu)^2/((x-xu)^2 + (y-yu)^2 + (z-zu)^2))^(1/2)))

