% This is the main script to be excuted for the first mechanism model
% (Pendulum with revolute joint) of Assignment 6 Problem 3

clear

% build the model
sys = simEngine3D("revJoint"); % revJoint.m is the model info file

t_start = 0;
t_end = 10;
t_step = 1e-3;

% Kinematics Analysis for 10 sec
states = sys.kinematic_analysis(t_start,t_step,t_end);

% separate r_i from q_i
r_i = zeros(3,length(states));
dr_i = zeros(3,length(states));
ddr_i = zeros(3,length(states));
for k = 1:length(states)
    r_i(:,k) = states{k}.q(1:3); % position of O'
    dr_i(:,k) = states{k}.dq(1:3); % velocity of O'
    ddr_i(:,k) = states{k}.ddq(1:3); % acceleration of O'
end

t = t_start:t_step:t_end;
% plot for point O'
figure 
subplot(3,1,1)
hold on
plot(t,r_i)
title('Position analysis of point O`')
xlabel('Time (sec)')
ylabel('Position in G-LF (m)')
legend('ri_X','ri_Y','ri_Z')
hold off

subplot(3,1,2)
hold on
plot(t,dr_i)
title('Velocity analysis of point O`')
xlabel('Time (sec)')
ylabel('Velocity in G-LF (m/s)')
legend('dri_X','dri_Y','dri_Z')
hold off

subplot(3,1,3)
hold on
plot(t,ddr_i)
title('Acceleration analysis of point O`')
xlabel('Time (sec)')
ylabel('Acceleration in G-LF (m/s^2)')
legend('ddri_X','ddri_Y','ddri_Z')
hold off

% plot for point Q
s_i_Q_bar = [-2, 0, 0]';
r_Q = zeros(3,length(states));
dr_Q = zeros(3,length(states));
ddr_Q = zeros(3,length(states));


