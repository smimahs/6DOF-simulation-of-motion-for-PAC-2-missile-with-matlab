clc;
clear all;
close all;

%% -------- Inputs --------
global w mu re J m_dot Thrust e1 e2 J1 J2 X_BC acceleration_g acceleration_nong l d i Ixx Iyy Izz m acceleration_body FA MA...
    Tempreture a_atmosisa P rho M0 m_fuel t_end_thrust i_end;

t_end_thrust = -1;
i_end = -1;
t = [0:0.1:350]; % simulation time
mu = 3.98603684e14;
re = 6378165.0; %m
J = 1.0823e-3;
e1 = 0;
e2 = 0;
J1 = 0;
J2 = 0;
X_BC = 1;

%% initial euler angle
phi = 0;
theta = 70;
psi = 90;

phi = deg2rad(phi);
theta = deg2rad(theta);
psi = deg2rad(psi);

%% initial position
lat = 24.7;
long = 46.7;
altitude = 620; 

R_I0 = lla2eci([lat,long,altitude],[2021 1 1 12 00 00]); % position in Inertia system
%R_I0 = -R_I0;
norm(R_I0);

% w
p = 0;
q = 0;
r = 0;
p = deg2rad(p);
q = deg2rad(q);
r = deg2rad(r);
w = [p,q,r];

%Thrust = 2000;
m_dot = 12.5; %kg/s
m0_total = 906; % kg
m_fuel = 420; % kg
M0 = m0_total - m_fuel;    
    
V0 = [0 0 0];
alpha = 0; % degree
alpha = deg2rad(alpha); % radian

l = 5.30; % length
d = 0.38608; % diameter

g0 = 9.81; 
i = 1; % counter

%% initial quaternion 

quaternion0 = angle2quat(psi,theta,phi);
quaternion0 = quatnormalize(quaternion0); % normalize quaternion

y0 = [R_I0 V0 w quaternion0];
[T,y] = ode23(@calculate,t,y0);

i_end;
t_end_thrust;
%% find quaternion magnetitude
%% find euler angles
for j=1:1:length(T)
    r_eci(j,:) = eci2lla(y(j,1:3),[2021 1 1 12 00 00]);
 
    V_norm(j) = norm(y(j,4:6));
    r_norm(j) = norm(y(j,1:3));
end

for j = 1:1:length(acceleration_g(:,4))
    acceleration_g(j,5) = norm(acceleration_g(j,1:3));
end
%% Show Acceleration
figure;
subplot(4,2,1)
plot(acceleration_g(:,4),acceleration_g(:,1),acceleration_g(:,4),acceleration_g(:,2),acceleration_g(:,4),acceleration_g(:,3));
legend({'ax','ay','az'},'Location','southwest');
grid on;
title('Acceleration G')
xlabel('t');
ylabel('a');

subplot(4,2,2)
plot(acceleration_g(:,4),acceleration_g(:,5));
grid on;
title('Acceleration G norm')
xlabel('t');
ylabel('a');

subplot(4,2,3)
plot(acceleration_nong(:,4),acceleration_nong(:,1),acceleration_nong(:,4),acceleration_nong(:,2),acceleration_nong(:,4),acceleration_nong(:,3));
legend({'ax','ay','az'},'Location','southwest');
grid on;
title('Acceleration non G')
xlabel('t');
ylabel('a');

subplot(4,2,4)
plot(acceleration_body(:,4),acceleration_body(:,1),acceleration_body(:,4),acceleration_body(:,2),acceleration_body(:,4),acceleration_body(:,3));
legend({'ax','ay','az'},'Location','southwest');
grid on;
title('Acceleration Body')
xlabel('t');
ylabel('a');

%% Mass
figure;
subplot(2,1,1)
plot(m(:,2),m(:,1));
grid on;
title('Mass')
xlabel('t');
ylabel('m');

%% Thrust
subplot(2,1,2)
plot(Thrust(:,2),Thrust(:,1));
grid on;
title('Thrust')
xlabel('t');
ylabel('T');

%% position in inertia
figure;
subplot(2,2,1)
plot(T,-y(:,1));
legend({'x'},'Location','southwest');
grid on;
title('position');
xlabel('t');
ylabel('rx');

subplot(2,2,2)
plot(T,-y(:,2));
legend({'y'},'Location','southwest');
grid on;
title('position');
xlabel('t');
ylabel('ry');

subplot(2,2,3)
plot(T,-y(:,3));
legend({'z'},'Location','southwest');
grid on;
title('position');
xlabel('t');
ylabel('rz');

subplot(2,2,4)
plot(T,r_norm);
legend({'r'},'Location','southwest');
grid on;
title('position norm');
xlabel('t');
ylabel('r');


%% Velocity
figure;
subplot(2,1,1)
plot(T,y(:,4),T,-y(:,5),T,y(:,6));
legend({'Vx','Vy','Vz'},'Location','southwest');
% grid on;
title('velocity');
xlabel('t');
ylabel('V');

subplot(2,1,2)
plot(T,V_norm);
legend({'V'},'Location','southwest');
grid on;
title('velocity norm');
xlabel('t');
ylabel('V');



%% position in eci system
figure;
subplot(2,2,1)
plot(T,r_eci(:,1));
legend({'lat'},'Location','southwest');
grid on;
title('lat');
xlabel('t');
ylabel('lat');

subplot(2,2,2)
plot(T,r_eci(:,2));
legend({'long'},'Location','southwest');
grid on;
title('long');
xlabel('t');
ylabel('long');

subplot(2,2,3)
plot(T,r_eci(:,3));
legend({'h'},'Location','southwest');
grid on;
title('h');
xlabel('t');
ylabel('h');


%% w
figure;
plot(T,rad2deg(y(:,7)),T,rad2deg(y(:,8)),T,rad2deg(y(:,9)));
legend({'wx','wy','wz'},'Location','southwest');
grid on;
title('w');
xlabel('t');
ylabel('w');

%% normalize quaternion
for j=1:length(y)
    y(j,10:13)=quatnormalize(y(j,10:13));
end

%% quaternion and euler angles
for j=1:1:length(T)
    eulZYX(j,:) = rad2deg(quat2eul(y(j,10:13)));
    mag(j) = quatnorm(y(j,10:13));
end

figure;
subplot(3,1,1);
plot(T,y(:,10),T,y(:,11),T,y(:,12),T,y(:,13));
legend({'a','b','c','d'},'Location','southwest');
grid on;
title('quaternion');
xlabel('t');
ylabel('q');

subplot(3,1,2);
plot(T,eulZYX);
legend({'Psi','Theta','Phi'},'Location','southwest');
grid on;
xlabel('t');
ylabel('Euler angle');

subplot(3,1,3);
plot(T,mag);

%% inertia
figure;
subplot(3,1,1);
plot(Ixx(:,2),Ixx(:,1));
legend({'Ixx'},'Location','southwest');
grid on;
title('Ixx')
xlabel('t');
ylabel('Ixx');

subplot(3,1,2);
plot(Iyy(:,2),Iyy(:,1));
legend({'Iyy'},'Location','southwest');
grid on;
title('Iyy')
xlabel('t');
ylabel('Iyy');
subplot(3,1,3);
plot(Izz(:,2),Izz(:,1));
legend({'Izz'},'Location','southwest');
grid on;
title('Izz')
xlabel('t');
ylabel('Izz');

%% Aerodynamic forces
figure;
subplot(2,1,1);
plot(FA(:,4),FA(:,1),FA(:,4),FA(:,2),FA(:,4),FA(:,3));
legend({'Fx','Fy','Fz'},'Location','southwest');
grid on;
title('FA')
xlabel('t');
ylabel('FA');

subplot(2,1,2);
plot(MA(:,4),MA(:,1),MA(:,4),MA(:,2),MA(:,4),MA(:,3));
legend({'Mx','My','Mz'},'Location','southwest');
grid on;
title('MA')
xlabel('t');
ylabel('MA');

%% atmospher
figure;
subplot(2,2,1);
plot(Tempreture(:,2),Tempreture(:,1));
grid on;
xlabel('H');
ylabel('T');
subplot(2,2,2);
plot(P(:,2),P(:,1));
grid on;
xlabel('H');
ylabel('P');
subplot(2,2,3);
plot(rho(:,2),rho(:,1));
grid on;
xlabel('H');
ylabel('rho');

subplot(2,2,4);
plot(a_atmosisa(:,2),a_atmosisa(:,1));
grid on;
xlabel('H');
ylabel('a');

% h_start = r_eci(end,3);
% v = norm(y(end,4:6));
% theta_angle = deg2rad(eulZYX(end,2));
% g = 9.81;
% 
% max_flight_time = (v * sin(theta_angle) + sqrt((v * sin(theta_angle))^2 + 2 * g * h_start)) / g;
% 
% hmax = h_start + v^2 * sin(theta_angle)^2 / (2 * g);
% renge_x = v * cos(theta_angle) * (v * sin(theta_angle) + sqrt(((v * sin(theta_angle))^2) + 2*g*h_start))/g;
% x = [1:1000:renge_x];
% H = h_start + x.*tan(theta_angle) - g*(x.^2)/(2*(v^2)*(cos(theta_angle)^2));
% vi = sqrt((v*cos(theta_angle))^2 - (2.*g.*(H - h_start)))
% figure;
% subplot(2,1,1);
% plot(x,H);
% legend({'x','h'},'Location','southwest');
% title(['max height is ' num2str(hmax/1000) ' km']);
% grid on;
% xlabel('x');
% ylabel('h');
% subplot(2,1,2);
% plot([1:max_flight_time],vi(1:fix(max_flight_time)));