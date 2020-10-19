function G = p2G(p)
% This function build the G matrix from Euler parameter p

e = p(2:4);
if size(e,2)>1
    e=e';
end

e_0 = p(1);

G = [-e, -cross_matrix(e)+e_0*eye(3)];

ened