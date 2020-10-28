function R = rot_z(angle)
% rotation matrix of z axis
R = [cos(angle) -sin(angle) 0;
     sin(angle) cos(angle) 0;
     0 0 1];
end