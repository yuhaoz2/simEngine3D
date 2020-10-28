% This is the main script to be excuted for the first mechanism model
% (Pendulum with revolute joint) of Assignment 8 Problem 1

clear

% build the model
sys = simEngine3D("revJoint"); % revJoint.m is the model info file

% Assign start time, ending time and time step
t_start = 0;
t_end = 10;
t_step = 1e-3;
BDF_order = 2;

% Perform the dynamic analysis
start = tic;
results = sys.dynamic_analysis(t_start,t_step,t_end,BDF_order);
run_time = toc(start);

disp(['The simulation time is ' num2str(run_time) ' seconds, and step size is ' num2str(t_step)]);

t = t_start:t_step:t_end;

%% Plots

% Plot the torque

torque =zeros(3,length(t));

for i = 1: length(t)
    for j = 1:5
        torque(:,i) = torque(:,i) + results{i}.torques{j};
    end
end

figure
plot(t,torque);
title('Plot of torque as a function of time');
xlabel('Time (sec)');
ylabel('Torque (Nm)');
legend('T_z','T_y','T_x');

% Plot coordiantes of point O'

r = zeros(3,length(results));
for k = 1:length(results)
    r(:,k) = results{k}.r; % position of O'
end

figure 
subplot(3,1,1)
plot(t,r(1,:));
title('x coordinate function plot of point O`');
xlabel('Time (sec)');
ylabel('x coordinate (m)');

subplot(3,1,2)
plot(t,r(2,:));
title('y coordinate function plot of point O`');
xlabel('Time (sec)');
ylabel('y coordinate (m)');

subplot(3,1,3)
plot(t,r(3,:));
title('z coordinate function plot of point O`');
xlabel('Time (sec)');
ylabel('z coordinate (m)');

% Plot the 2-norm of the violation of the velocity constraint equations

vio = zeros(1,length(results));
for k = 1:length(results)
    vio(k) = norm(results{k}.violation_vel);
end

figure
plot(t,vio);
title('2-norm of the violation of the velocity constraint equations');
xlabel('Time (sec)');
ylabel('2-norm of violation');

% Plot the angular velocity of Body 1

p = zeros(4,length(results));
dp = zeros(4,length(results));
omega = zeros(4,length(results));

for k = 1:length(results)
    % compute orientation matrix
    p     = results{k}.p; 
    dp    = results{k}.dp; 

    omega(:,k) = 2*p2E(p)*dp;
end

figure 
subplot(3,1,1)
plot(t,omega(1,:));
title('Plot of angular velocity in G-RF x axis');
xlabel('Time (sec)');
ylabel('omega_x (rad/s)');

subplot(3,1,2)
plot(t,omega(2,:));
title('Plot of angular velocity in G-RF y axis');
xlabel('Time (sec)');
ylabel('omega_y (rad/s)');

subplot(3,1,3)
plot(t,omega(3,:));
title('Plot of angular velocity in G-RF z axis');
xlabel('Time (sec)');
ylabel('omega_z (rad/s)');
