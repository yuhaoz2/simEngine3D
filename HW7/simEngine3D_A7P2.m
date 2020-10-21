% This is the main script to be excuted for Assignment 7 Problem 2
% Use Backward Euler to find an approximation of the exact solution of this IVP.

clear

% Define the time duration
t_start = 0;
t_end = 20;

% step size
h = 0.001;

counts = length(t_start:h:t_end);

%  Initial value
x_0 = 0;
y_0 = 2;

% Initialize
x = zeros(counts,1);
y = zeros(counts,1);
x(1) = x_0;
y(1) = y_0;
eps = 1e-6; % tolerance
itCountMax = 50; % max iteration

for i=2:counts
    
    % start point
    x_i = x(i-1);
    y_i = y(i-1);
    
    for j = 1:itCountMax
        % Step 1: compute correction
        
        % jacobian
        J = [1+h+4*h*y_i*(1-x_i^2)/(1+x_i^2)^2 4*h*x_i/(1+x_i^2);...
            -h+h*(1-x_i^2)/(1+x_i^2)^2 1+h*x_i/(1+x_i^2)];
        % g(x,y)
        g = [x_i*(1+h)+4*h*x_i*y_i/(1+x_i^2)-x(i-1)-h;
            -h*x_i+y_i+h*x_i*y_i/(1+x_i^2)-y(i-1)];
        
        delta = J\-g;
        
        % Step 2: apply correction
        x_i = x_i + delta(1);
        y_i = y_i + delta(2);
        
        % Step 3: compute norm residual
        res = norm(delta);
        
        % Step 4
        if res < eps
            break;
        end
    end
    
    if j >= itCountMax
        error('Maximum number of iterations reached');
    end
    x(i) = x_i;
    y(i) = y_i;
    
end

% Plots
t=t_start:h:t_end;
figure
subplot(2,1,1)
plot(t,x);
title('Plot of x(t)');
xlabel('t');
ylabel('x');

subplot(2,1,2)
plot(t,y);
title('Plot of y(t)');
xlabel('t');
ylabel('y');


