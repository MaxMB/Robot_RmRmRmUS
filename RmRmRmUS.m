close, clear, clc;
%%% SYSTEM CONSTANTS
l1 = 0.2;
%lG1 = l1/2;
l2 = 0.2;
lG2 = l2/2;
l3 = 1;
lG3 = l3/2;
l4 = 1;
lB = 1;
%lGB = lB/2;
x0S = 0;
y0S = 2.5;
z0S = 0.5;
J1 = 0.1;
J2 = 0.1;
J3 = 0.1;
%JBx = 0.1;
%JBz = 0.1;
%m1 = 1;
m2 = 1;
m3 = 1;
%mB = 1;
g = 9.8;


%%% INVERSE KINEMATICS
% Initial and final conditions
qB0 = [deg2rad(185); deg2rad(10)]; % qB(0) = [q4(0); q5(0)]
qBf = [deg2rad(175); deg2rad(20)]; % qB(tf) = [q4(tf); q5(tf)]
% Polynomial initial and final values
theta0 = 0;
thetaf = norm(qBf-qB0);
dtheta0 = 0;
dthetaf = 0;
ddtheta0 = 0;
ddthetaf = 0;
ndv = (qBf-qB0) / thetaf; % normalized direction vector
% Quintic polynomial
tf_ik = 1;
a0 = theta0;
a1 = dtheta0;
a2 = ddtheta0/2;
a3 = ( 20*(thetaf - theta0) ...
      - (8*dthetaf + 12*dtheta0)*tf_ik ...
      - (3*ddtheta0 - ddthetaf)*tf_ik^2 ) / (2*tf_ik^3);
a4 = ( -30*(thetaf - theta0) ...
      + (14*dthetaf + 16*dtheta0)*tf_ik ...
      + (3*ddtheta0 - 2*ddthetaf)*tf_ik^2 ) / (2*tf_ik^4);
a5 = ( 12*(thetaf - theta0) ...
      - (6*dthetaf + 6*dtheta0)*tf_ik ...
      - (ddtheta0 - ddthetaf)*tf_ik^2 ) / (2*tf_ik^5);
% Cubic polynomial
%a0 = theta0;
%a1 = 0;
%a2 = 3*(thetaf-theta0)/tf_ik^2;
%a3 = -2*(thetaf-theta0)/tf_ik^3;


%%% INITIAL CONDITIONS
% Bar direct kinematics
pxT_0 = - cos(qB0(2)) * sin(qB0(1)) * lB + x0S;
pyT_0 = cos(qB0(2)) * cos(qB0(1)) * lB + y0S;
pzT_0 = sin(qB0(2)) * lB + z0S;
% Robotic arm inverse kinematics
dzT_0 = pzT_0 - l1 - l2;
Rxy_0 = sqrt(pxT_0^2 + pyT_0^2);
q1_0 = asin(-pxT_0/Rxy_0);
q3_0 = -acos((Rxy_0^2+dzT_0^2-l3^2-l4^2)/(2*l3*l4));
q2_0 = -asin(((l3+l4*cos(q3_0))*Rxy_0+l4*sin(q3_0)*dzT_0)/...
    (l3^2+2*l3*l4*cos(q3_0)+l4^2));
% Initial
qI0 = [q1_0; q2_0];
qD0 = [q3_0; qB0];
dqI0 = [0; 0];
T0 = [0; 0; 0];


%%% CONTROLLER
t_stl = tf_ik;
Kp = [1; 1; 1] * (4/t_stl)^2;
Kd = [1; 1; 1] * (8/t_stl);


%%% SIMULATION
tf = tf_ik + 0.5;
t_step = 1e-4;
sim('RmRmRmUS_Simulink.slx');
% Generalized coordinates
q1 = qI(:,1);
q2 = qI(:,2);
q3 = qD(:,1);
q4 = qD(:,2);
q5 = qD(:,3);
% Generalized velocities
dq1 = dqI(:,1);
dq2 = dqI(:,2);
dq3 = dqD(:,1);
dq4 = dqD(:,2);
dq5 = dqD(:,3);
% Torques
T1 = To(:,1);
T2 = To(:,2);
T3 = To(:,3);
% Tool position by the bar
xT_B = pT(:,1);
yT_B = pT(:,2);
zT_B = pT(:,3);
% Tool position by the robotic arm
xT_R = sin(q1) .* (l3*sin(q2) + l4*sin(q2+q3));
yT_R = - cos(q1) .* (l3*sin(q2) + l4*sin(q2+q3));
zT_R = l1 + l2 + l3*cos(q2) + l4*cos(q2+q3);
% Reference generalized coordinates
qR1 = qR(:,1);
qR2 = qR(:,2);
qR3 = qR(:,3);
% Reference generalized velocities
dqR1 = dqR(:,1);
dqR2 = dqR(:,2);
dqR3 = dqR(:,3);


%%% Plot
figure(1), set(gcf,'color','w');
subplot(511), plot(t,rad2deg(qR1),'r'), hold on;
    plot(t,rad2deg(q1),'--k'), hold off;
    grid on, ylabel('q_1'), legend('Reference','Robot');
    title('Joint Position [º]');
subplot(512), plot(t,rad2deg(qR2),'r'), hold on;
    plot(t,rad2deg(q2),'--k'), hold off, grid on, ylabel('q_2');
subplot(513), plot(t,rad2deg(qR3),'r'), hold on;
    plot(t,rad2deg(q3),'--k'), hold off, grid on, ylabel('q_3');
subplot(514), plot(t,rad2deg(q4),'k'), grid on, ylabel('q_4');
    legend('Bar');
subplot(515), plot(t,rad2deg(q5),'k'), grid on, ylabel('q_5');
    xlabel('Time [s]');

figure(2), set(gcf,'color','w');
subplot(511), plot(t,rad2deg(dqR1),'r'), hold on;
    plot(t,rad2deg(dq1),'--k'), hold off;
    grid on, ylabel('dq_1/dt'), legend('Reference','Robot');
    title('Joint Velocity [º/s]');
subplot(512), plot(t,rad2deg(dqR2),'r'), hold on;
    plot(t,rad2deg(dq2),'--k'), hold off, grid on, ylabel('dq_2/dt');
subplot(513), plot(t,rad2deg(dqR3),'r'), hold on;
    plot(t,rad2deg(dq3),'--k'), hold off, grid on, ylabel('dq_3/dt');
subplot(514), plot(t,rad2deg(dq4),'k'), grid on, ylabel('dq_4/dt');
    legend('Bar');
subplot(515), plot(t,rad2deg(dq5),'k'), grid on, ylabel('dq_5/dt');
    xlabel('Time [s]');

figure(3), set(gcf,'color','w');
subplot(311), plot(t,T1,'k'), grid on, ylabel('T_1');
    title('Joint Torque [N*m]');
subplot(312), plot(t,T2,'k'), grid on, ylabel('T_2');
subplot(313), plot(t,T3,'k'), grid on, ylabel('T_3'), xlabel('Time [s]');

figure(4), set(gcf,'color','w');
plot3(xT_B, yT_B, zT_B, 'r'), hold on;
plot3(xT_B(1), yT_B(1), zT_B(1), 'sb');
plot3(xT_R, yT_R, zT_R, '--k');
plot3(xT_R(1), yT_R(1), zT_R(1), '*g');
plot3(x0S, y0S, z0S, 'om');
plot3(0, 0, 0,'oc'), hold off;
xlabel('x [m]'), ylabel('y [m]'), zlabel('z [m]'), grid on;
legend('Bar Trajectory', 'Init Bar', 'Robot Trajectory', 'Init Robot',...
    'Origin Bar', 'Origin Robot');