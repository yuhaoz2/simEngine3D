function A = p2A(p)
% This function takes a Euler parameter p and convert it into the orientation matrix A

e = p(2:4);
if size(e,2)>1
    e=e';
end

A = (2*p(1)^2-1)*eye(3) + 2*(e*e'+p(1)*cross_matrix(e));
end