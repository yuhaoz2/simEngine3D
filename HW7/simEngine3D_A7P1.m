% This is the main script to be excuted for the first mechanism model
% (Pendulum with revolute joint) of Assignment 7 Problem 1

clear

% build the model
sys = simEngine3D("revJoint"); % revJoint.m is the model info file

% Assign start time, ending time and time step
t_start = 0;
t_end = 10;
t_step = 1e-3;

% Perform the inverse dynamic analysis
torques = sys.inverse_dynamics(t_start,t_step,t_end);

t = t_start:t_step:t_end;

torque =zeros(3,length(t));

% Find the torque on the driving constraint
for i = 1: length(t)
    torque(:,i) = torques{i}{6,1};
end

% Plot the torque
plot(t,torque);
title('Torque applied to pendulum');
xlabel('Time (sec)');
ylabel('Torque (Nm)');
legend('T_z','T_y','T_x');