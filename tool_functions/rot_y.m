function R = rot_y(angle)
% rotation matrix of y axis
R = [cos(angle) 0 sin(angle);
     0 1 0;
     -sin(angle) 0 cos(angle)];
end