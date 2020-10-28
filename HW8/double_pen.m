% This script define the double pendulum model for HW8 Problem 2

% Define number of bodies
nb = 2;
g = [0;0;-9.81];  


% body 1
body(1).id = 1;
body(1).r_0 = [0, 2, 0]';
body(1).p_0 = A2p(rot_y(pi/2)*rot_z(pi/2));
body(1).dr_0 = [0, 0, 0]';
body(1).dp_0 = [0, 0, 0, 0]';
body(1).ddr_0 = [0, 0, 0]';
body(1).ddp_0 = [0, 0, 0, 0]';

body(1).mass = 7800*4*0.05*0.05;
body(1).inertia = [body(1).mass/12*2*0.05^2 0 0;...
                   0 body(1).mass/12*(4^2+0.05^2) 0;...
                   0 0 body(1).mass/12*(4^2+0.05^2)];
             
body(1).force = body(1).mass*g;
body(1).torque = [0 0 0]';

% body 2
body(2).id = 2;
body(2).r_0 = [0, 4, -1]';
body(2).p_0 = A2p(rot_y(pi/2));
body(2).dr_0 = [0, 0, 0]';
body(2).dp_0 = [0, 0, 0, 0]';
body(2).ddr_0 = [0, 0, 0]';
body(2).ddp_0 = [0, 0, 0, 0]';

body(2).mass = 7800*2*0.05*0.05;
body(2).inertia = [body(2).mass/12*2*0.05^2 0 0;...
                   0 body(2).mass/12*(2^2+0.05^2) 0;...
                   0 0 body(2).mass/12*(2^2+0.05^2)];
             
body(2).force = body(2).mass*g;
body(2).torque = [0 0 0]';

% Define number of constaints
nc = 10;

% constraint between body 1 and ground

% 3 CD constraints
% constraint 1
constraints(1).type = "CD";
constraints(1).body_i = 1;
constraints(1).body_j = 0;
constraints(1).c = [1, 0, 0]';
constraints(1).s_i_P_bar = [-2, 0, 0]';
constraints(1).s_j_Q_bar =[0, 0, 0]';
constraints(1).f = @(t)(0*t);
constraints(1).df = @(t)(0*t);
constraints(1).ddf = @(t)(0*t);

% constraint 2
constraints(2).type = "CD";
constraints(2).body_i = 1;
constraints(2).body_j = 0;
constraints(2).c = [0, 1, 0]';
constraints(2).s_i_P_bar = [-2, 0, 0]';
constraints(2).s_j_Q_bar =[0, 0, 0]';
constraints(2).f = @(t)(0*t);
constraints(2).df = @(t)(0*t);
constraints(2).ddf = @(t)(0*t);

% constraint 3
constraints(3).type = "CD";
constraints(3).body_i = 1;
constraints(3).body_j = 0;
constraints(3).c = [0, 0, 1]';
constraints(3).s_i_P_bar = [-2, 0, 0]';
constraints(3).s_j_Q_bar =[0, 0, 0]';
constraints(3).f = @(t)(0*t);
constraints(3).df = @(t)(0*t);
constraints(3).ddf = @(t)(0*t);

% 2 DP1 constraints
% constraint 4
constraints(4).type = "DP1";
constraints(4).body_i = 1;
constraints(4).body_j = 0;
constraints(4).a_i_bar = [1, 0, 0]';
constraints(4).a_j_bar =[1, 0, 0]';
constraints(4).f = @(t)(0*t);
constraints(4).df = @(t)(0*t);
constraints(4).ddf = @(t)(0*t);

% constraint 5
constraints(5).type = "DP1";
constraints(5).body_i = 1;
constraints(5).body_j = 0;
constraints(5).a_i_bar = [0, 1, 0]';
constraints(5).a_j_bar =[1, 0, 0]';
constraints(5).f = @(t)(0*t);
constraints(5).df = @(t)(0*t);
constraints(5).ddf = @(t)(0*t);

% constraint between body 1 and body 2

% 3 CD constraints
% constraint 6
constraints(6).type = "CD";
constraints(6).body_i = 1;
constraints(6).body_j = 2;
constraints(6).c = [1, 0, 0]';
constraints(6).s_i_P_bar = [2, 0, 0]';
constraints(6).s_j_Q_bar =[-1, 0, 0]';
constraints(6).f = @(t)(0*t);
constraints(6).df = @(t)(0*t);
constraints(6).ddf = @(t)(0*t);

% constraint 7
constraints(7).type = "CD";
constraints(7).body_i = 1;
constraints(7).body_j = 2;
constraints(7).c = [0, 1, 0]';
constraints(7).s_i_P_bar = [2, 0, 0]';
constraints(7).s_j_Q_bar =[-1, 0, 0]';
constraints(7).f = @(t)(0*t);
constraints(7).df = @(t)(0*t);
constraints(7).ddf = @(t)(0*t);

% constraint 8
constraints(8).type = "CD";
constraints(8).body_i = 1;
constraints(8).body_j = 2;
constraints(8).c = [0, 0, 1]';
constraints(8).s_i_P_bar = [2, 0, 0]';
constraints(8).s_j_Q_bar =[-1, 0, 0]';
constraints(8).f = @(t)(0*t);
constraints(8).df = @(t)(0*t);
constraints(8).ddf = @(t)(0*t);

% 2 DP1 constraints
% constraint 9
constraints(9).type = "DP1";
constraints(9).body_i = 1;
constraints(9).body_j = 2;
constraints(9).a_i_bar = [0, 0, 1]';
constraints(9).a_j_bar =[1, 0, 0]';
constraints(9).f = @(t)(0*t);
constraints(9).df = @(t)(0*t);
constraints(9).ddf = @(t)(0*t);

% constraint 10
constraints(10).type = "DP1";
constraints(10).body_i = 1;
constraints(10).body_j = 2;
constraints(10).a_i_bar = [0, 0, 1]';
constraints(10).a_j_bar =[0, 1, 0]';
constraints(10).f = @(t)(0*t);
constraints(10).df = @(t)(0*t);
constraints(10).ddf = @(t)(0*t);