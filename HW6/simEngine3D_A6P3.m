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

for k = 1:length(states)
    % compute orientation matrix
    p_i     = states{k}.q(4:7); 
    dp_i    = states{k}.dq(4:7); 
    ddp_i   = states{k}.ddq(4:7); 
    A_i = p2A(p_i);
    B = pa2B(p_i,s_i_Q_bar);
    r_Q(:,k) = r_i(:,k) + A_i*s_i_Q_bar;
    dr_Q(:,k) = dr_i(:,k) + B*dp_i;
    ddr_Q(:,k) = ddr_i(:,k) + B*ddp_i + pa2B(dp_i,s_i_Q_bar)*dp_i;
end

figure
subplot(3,1,1)
hold on
plot(t,r_Q)
title('Position analysis of point Q')
xlabel('Time (sec)')
ylabel('Position in G-LF (m)')
legend('rQ_X','rQ_Y','rQ_Z')
hold off

subplot(3,1,2)
hold on
plot(t,dr_Q)
title('Velocity analysis of point Q')
xlabel('Time (sec)')
ylabel('Velocity in G-LF (m/s)')
legend('drQ_X','drQ_Y','drQ_Z')
hold off

subplot(3,1,3)
hold on
plot(t,ddr_Q)
title('Acceleration analysis of point Q')
xlabel('Time (sec)')
ylabel('Acceleration in G-LF (m/s^2)')
legend('ddrQ_X','ddrQ_Y','ddrQ_Z')
hold off

