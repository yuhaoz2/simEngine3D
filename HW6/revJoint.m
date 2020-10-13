% body 1
i = 1;
r_i_0 = [0, sqrt(2), -sqrt(2)]';
p_i_0 = [0.6533, 0.2706, 0.6533, 0.2706]';
q_i_0 = [r_i_0; p_i_0];
dr_i_0 = [0, 0, 0]';
dp_i_0 = [0, 0, 0, 0]';
dq_i_0 = [dr_i_0; dp_i_0];

% ground
j = 0;
q_j_0 = [];
dq_j_0 = [];

syms f(t);

% 3 CD constraints
% constraint 1
c_1 = [1, 0, 0]';
s_i_P_bar_1 = [-2, 0, 0]';
s_j_Q_bar_1 =[0, 0, 0]';
f_1 = 0;
df_1 = 0;
ddf_1 = 0;

% constraint 2
c_2 = [0, 1, 0]';
s_i_P_bar_2 = [-2, 0, 0]';
s_j_Q_bar_2 =[0, 0, 0]';
f_2 = 0;
df_2 = 0;
ddf_2 = 0;

% constraint 3
c_3 = [0, 0, 1]';
s_i_P_bar_3 = [-2, 0, 0]';
s_j_Q_bar_3 =[0, 0, 0]';
f_3 = 0;
df_3 = 0;
ddf_3 = 0;

% 2 DP1 constraints
% constraint 4
a_i_bar_4 = [1, 0, 0]';
a_j_bar_4 =[1, 0, 0]';
f_4 = 0;
df_4 = 0;
ddf_4 = 0;

% constraint 5
a_i_bar_5 = [0, 1, 0]';
a_j_bar_5 =[1, 0, 0]';
f_5 = 0;
df_5 = 0;
ddf_5 = 0;

% driving constraint
a_i_bar_6 = [1, 0, 0]';
a_j_bar_6 =[0, 0, -1]';
f(t) = cos(pi/4*cos(2*t));
df = diff(f,t);
ddf = diff(df,t);