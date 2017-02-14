function ceq = p3_nonlinear_couple( A, b, c, At, bt, ct, isPim2)
% function ceq = p3_nonlinear_couple( A, b, c, At, bt, ct)
% Purpose: ARK third-order nonlinear coupled conditions

A; At;

ceq(1) = b'*(c.*ct) - 1/3;
ceq(2) = b'*(ct.*c) - 1/3;
ceq(3) = b'*(ct.*ct) - 1/3;
if isPim2
    ceq(4) = 0; %bt'*(c.*c) - 1/3;
    ceq(5) = 0; %bt'*(c.*ct) - 1/3;
    ceq(6) = 0; %bt'*(ct.*c) - 1/3;
else
    ceq(4) = bt'*(c.*c) - 1/3;
    ceq(5) = bt'*(c.*ct) - 1/3;
    ceq(6) = bt'*(ct.*c) - 1/3;
end

end
