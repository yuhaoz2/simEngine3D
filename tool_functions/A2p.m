function p = A2p(A)
% convert orentation matrix A to euler parameter p
e0 = sqrt((trace(A)+1)/4);

if e0 ~= 0
    e1 = (A(3,2)-A(2,3))/(4*e0);
    e2 = (A(1,3)-A(3,1))/(4*e0);
    e3 = (A(2,1)-A(1,2))/(4*e0);
else
    e1 = sqrt((2*A(1,1)-trace(A)+1)/4);
    if e1 == 0
        e2 = sqrt((2*A(2,2)-trace(A)+1)/4);
        if e2 == 0
            e3 = 1;
        else
            e3 = (A(3,2)+A(2,3))/(4*e2);
        end
    else
        e2 = (A(2,1)+A(1,2))/(4*e1);;
        e3 = (A(1,3)+A(3,1))/(4*e1);
    end
end

p = [e0;e1;e2;e3];

end