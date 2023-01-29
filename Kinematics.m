close, clear, clc;
%%% Initial and final conditions
qB0 = [deg2rad(185); deg2rad(10)]; % qB(0) = [q4(0); q5(0)]
qBf = [deg2rad(175); deg2rad(20)]; % qB(tf) = [q4(tf); q5(tf)]
% Polynomial initial and final values
theta0 = 0;
thetaf = norm(qBf-qB0);
dtheta0 = 0;
dthetaf = 0;
ddtheta0 = 0;
ddthetaf = 0;
% Time
t_step = 1e-4;
tf = 1;
t = 0:t_step:tf;
% Quintic polynomial
a0 = theta0;
a1 = dtheta0;
a2 = ddtheta0/2;
a3 = ( 20*(thetaf - theta0) ...
      - (8*dthetaf + 12*dtheta0)*tf ...
      - (3*ddtheta0 - ddthetaf)*tf^2 ) / (2*tf^3);
a4 = ( -30*(thetaf - theta0) ...
      + (14*dthetaf + 16*dtheta0)*tf ...
      + (3*ddtheta0 - 2*ddthetaf)*tf^2 ) / (2*tf^4);
a5 = ( 12*(thetaf - theta0) ...
      - (6*dthetaf + 6*dtheta0)*tf ...
      - (ddtheta0 - ddthetaf)*tf^2 ) / (2*tf^5);
theta = a0 + a1*t + a2*t.^2 + a3*t.^3 + a4*t.^4 + a5*t.^5;
dtheta = a1 + 2*a2*t + 3*a3*t.^2 + 4*a4*t.^3 + 5*a5*t.^4;
ddtheta = 2*a2 + 6*a3*t + 12*a4*t.^2 + 20*a5*t.^3;
% Cubic polynomial
% a0 = theta0;
% a1 = 0;
% a2 = 3*(thetaf-theta0)/tf^2;
% a3 = -2*(thetaf-theta0)/tf^3;
% theta = a0 + a1*t + a2*t.^2 + a3*t.^3;
% dtheta = a1 + 2*a2*t + 3*a3*t.^2;
% ddtheta = 2*a2 + 6*a3*t;
% Bar generalized coordinates q4 and q5
ndv = (qBf-qB0) / norm(qBf-qB0); % normalized direction vector
n = length(t);
qB = zeros(2,n);
dqB = zeros(2,n);
ddqB = zeros(2,n);
for i = 1:n
    qB(:,i) = qB0 + theta(i) * ndv; % Bar generalized coordinates
    dqB(:,i) = dtheta(i) * ndv; % Bar generalized velocities
    ddqB(:,i) = ddtheta(i) * ndv; % Bar generalized accelerations
end
% Bar generalized coordinates
q4 = qB(1,:);
q5 = qB(2,:);
% Bar generalized velocities
dq4 = dqB(1,:);
dq5 = dqB(2,:);
% Bar generalized accelerations
ddq4 = ddqB(1,:);
ddq5 = ddqB(2,:);



%%% Bar Direct Kinematics
lB = 1;
x0S = 0;
y0S = 2.5;
z0S = 0.5;
% Tool position
xT = - cos(q5) .* sin(q4) * lB + x0S;
yT = cos(q5) .* cos(q4) * lB + y0S;
zT = sin(q5) * lB + z0S;
% Tool velocity and acceleration
vT = zeros(3,n);
aT = zeros(3,n);
for i = 1:n
    JqB = mtx_JqB(q4(i), q5(i), lB);
    dJqB = mtx_dJqB(q4(i), q5(i), dq4(i), dq5(i), lB);
    vT(:,i) = JqB * dqB(:,i);
    aT(:,i) = dJqB * dqB(:,i) + JqB * ddqB(:,i);
end



%%% Robotic Arm Inverse Kinematics
l1 = 0.2;
l2 = 0.2;
l3 = 1;
l4 = 1;
dzT = zT - l1 - l2;
Rxy = sqrt(xT.^2 + yT.^2);
% Robotic arm generalized coordinates
q1 = asin(- xT ./ Rxy);
q3 = - acos((Rxy.^2 + dzT.^2 - l3^2 - l4^2) / (2*l3*l4));
q2 = - asin(((l3+l4*cos(q3)).*Rxy + l4*sin(q3).*dzT) ./ ...
    (l3^2 + 2*l3*l4*cos(q3) + l4^2));
% Robotic arm generalized velocities and accelerations
dqR = zeros(3,n);
ddqR = zeros(3,n);
for i = 1:n
    iJqR = mtx_invJqR(q1(i), q2(i), q3(i), l3, l4);
    dqR(:,i) = iJqR * vT(:,i);
    dJqR = mtx_dJqR(q1(i), q2(i), q3(i), dqR(1,i), dqR(2,i), dqR(3,i), l3, l4);
    ddqR(:,i) = iJqR * (aT(:,i) - dJqR * dqR(:,i));
end
% Robotic arm generalized velocities
dq1 = dqR(1,:);
dq2 = dqR(2,:);
dq3 = dqR(3,:);
% Robotic arm generalized accelerations
ddq1 = ddqR(1,:);
ddq2 = ddqR(2,:);
ddq3 = ddqR(3,:);



%%% OUTPUT
figure(1), set(gcf,'color','w');
subplot(511), plot(t,rad2deg(q1),'k'), grid on, ylabel('q_1 [º]');
    title('Inverse Kinematics 5º Polynomial - Generalized Coordinates');
subplot(512), plot(t,rad2deg(q2),'k'), grid on, ylabel('q_2 [º]');
subplot(513), plot(t,rad2deg(q3),'k'), grid on, ylabel('q_3 [º]');
subplot(514), plot(t,rad2deg(q4),'k'), grid on, ylabel('q_4 [º]');
subplot(515), plot(t,rad2deg(q5),'k'), grid on, ylabel('q_5 [º]');
    xlabel('Time [s]');

figure(2), set(gcf,'color','w');
subplot(511), plot(t,rad2deg(dq1),'k'), grid on, ylabel('dq_1/dt [º/s]');
    title('Inverse Kinematics 5º Polynomial - Generalized Velocities');
subplot(512), plot(t,rad2deg(dq2),'k'), grid on, ylabel('dq_2/dt [º/s]');
subplot(513), plot(t,rad2deg(dq3),'k'), grid on, ylabel('dq_3/dt [º/s]');
subplot(514), plot(t,rad2deg(dq4),'k'), grid on, ylabel('dq_4/dt [º/s]');
subplot(515), plot(t,rad2deg(dq5),'k'), grid on, ylabel('dq_5/dt [º/s]');
    xlabel('Time [s]');

figure(3), set(gcf,'color','w');
subplot(511), plot(t,rad2deg(ddq1),'k'), grid on;
    ylabel('d^2q_1/dt^2 [º/s^2]');
    title('Inverse Kinematics 5º Polynomial - Generalized Accelerations');
subplot(512), plot(t,rad2deg(ddq2),'k'), grid on;
    ylabel('d^2q_2/dt^2 [º/s^2]');
subplot(513), plot(t,rad2deg(ddq3),'k'), grid on;
    ylabel('d^2q_3/dt^2 [º/s^2]');
subplot(514), plot(t,rad2deg(ddq4),'k'), grid on;
    ylabel('d^2q_4/dt^2 [º/s^2]');
subplot(515), plot(t,rad2deg(ddq5),'k'), grid on;
    ylabel('d^2q_5/dt^2 [º/s^2]');
    xlabel('Time [s]');



%{
xT2 = sin(q1) .* (l3*sin(q2) + l4*sin(q2+q3));
yT2 = - cos(q1) .* (l3*sin(q2) + l4*sin(q2+q3));
zT2 = l1 + l2 + l3*cos(q2) + l4*cos(q2+q3);

xr2_0 = sin(q1(1)) * l3*sin(q2(1));
yr2_0 = - cos(q1(1)) * l3*sin(q2(1));
zr2_0 = l1 + l2 + l3*cos(q2(1));

xr3_0 = xr2_0 + sin(q1(1)) * l4*sin(q2(1)+q3(1));
yr3_0 = yr2_0 - cos(q1(1)) * l4*sin(q2(1)+q3(1));
zr3_0 = zr2_0 + l4*cos(q2(1)+q3(1));

xrT_0 = [0, xr2_0, xr3_0];
yrT_0 = [0, yr2_0, yr3_0];
zrT_0 = [0, zr2_0, zr3_0];

figure(1), set(gcf,'color','w');
plot3(xT, yT, zT, 'r'), hold on;
plot3(xT(1), yT(1), zT(1), 'sb');
plot3(xT2, yT2, zT2, '--k');
plot3(xT2(1), yT2(1), zT2(1), '*g');
plot3(x0S, y0S, z0S, 'om');
plot3(0, 0, 0,'oc');
plot3(xrT_0, yrT_0, zrT_0, 'k'), hold off;
xlabel('x [m]'), ylabel('y [m]'), zlabel('z [m]'), grid on;
legend('Bar Trajectory', 'Init Bar', 'Robot Trajectory', 'Init Robot',...
    'Origin Bar', 'Origin Robot', 'Robotic Arms');
%}