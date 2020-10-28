function R = rot_x(angle)
% rotation matrix of x axis
R = [1 0 0;
     0 cos(angle) -sin(angle);
     0 sin(angle) cos(angle)];
end