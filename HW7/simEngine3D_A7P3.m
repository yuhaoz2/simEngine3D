% This is the main script to be excuted for Assignment 7 Problem 3
% convergence analysis of an IVP.

clear


% Define the time duration
t_start = 1;
t_end = 10;

% step size
h = [0.001 0.002 0.004 0.005 0.008 0.009 0.01 0.015 0.02];

%  Initial value
y_0 = 1;

% ground trueth value at t=10
y_10 = 1/10+tan(1/10+pi-1)/10^2;

% First use Backward Euler

% record Backward Euler solution error
error_BE = zeros(length(h),1);

for k=1:length(h)
    
    % Initialize
    counts = length(t_start:h(k):t_end);
    t=t_start:h(k):t_end;
    y = zeros(counts,1);
    y(1) = y_0;
    eps = 1e-8; % tolerance
    itCountMax = 50; % max iteration
    
    for i=2:counts
        
        % start point
        y_i = y(i-1);
        
        for j = 1:itCountMax
            % Step 1: compute correction
            
            % jacobian
            J = 1 + 2*h(k)*y_i;
            % g(y)
            g = y_i-y(i-1)-h(k)*(-y_i^2-1/t(i)^4);
            
            delta = -g/J;
            
            % Step 2: apply correction
            y_i = y_i + delta;
            
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
        y(i) = y_i;
    end
    error_BE(k) = abs(y(end)-y_10);
end


% Secondly use 4th order BDF

% Define the coefficients
beta = 12/25;
a1 = -48/25;
a2 = 36/25;
a3 = -16/25;
a4 = 3/25;

% record SBDF Euler solution error
error_BDF = zeros(length(h),1);

for k=1:length(h)
    
    % Initialize
    counts = length(t_start:h(k):t_end);
    t=t_start:h(k):t_end;
    y = zeros(counts,1);
    y(1) = y_0;
    y(2) = 1/t(2)+tan(1/t(2)+pi-1)/t(2)^2;
    y(3) = 1/t(3)+tan(1/t(3)+pi-1)/t(3)^2;
    y(4) = 1/t(4)+tan(1/t(4)+pi-1)/t(4)^2;
    
    eps = 1e-10; % tolerance
    itCountMax = 50; % max iteration
    
    for i=5:counts
        
        % start point
        y_i = y(i-1);
        
        for j = 1:itCountMax
            % Step 1: compute correction
            
            % jacobian
            J = 1 + 2*beta*h(k)*y_i;
            % g(y)
            g = y_i+a1*y(i-1)+a2*y(i-2)+a3*y(i-3)+a4*y(i-4)-beta*h(k)*(-y_i^2-1/t(i)^4);
            delta = -g/J;
            
            % Step 2: apply correction
            y_i = y_i + delta;
            
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
        y(i) = y_i;
    end
    error_BDF(k) = abs(y(end)-y_10);
end

% Plots

figure
subplot(2,1,1)
plot(h,error_BE);
title('Backward Euler convergence plot');
ylabel('error (|y(t_n)-y_n|)');
xlabel('h');
subplot(2,1,2)
plot(h,error_BDF);
title('4th order BDF convergence plot plot');
ylabel('error (|y(t_n)-y_n|)');
xlabel('h');

figure
plot(log(h),log(error_BE),'r',log(h),log(error_BDF),'b');
title('Convergence plot (log)');
ylabel('log(error)');
xlabel('log(h)');