function ceq = p4_nonlinear_couple( A, b, c, At, bt, ct, isPim2)
% function ceq = p4_nonlinear_couple( A, b, c, At, bt, ct)
% Purpose: ARK Fourth-order nonlinear coupled conditions
% the line numbers corresponds to the the line number of the trees from the
% imex-order condition tree notes done by Sidafa Conde

% The conditions commented out are either:
%   - the non-coupled explicit nonlinear conditions
%   - the non-coupled implicit nonlinear conditions
%   - the coupled and non-coupled linear conditions
%       - and this is already considered in the general_coupled_linear

C = diag(c); Ct = diag(ct);

% line 1
ceq(1) = b'*(c.*c.*c) - 1/4;
ceq(2) = b'*(c.*c.*ct) - 1/4;
ceq(3) = b'*C*A*c - 1/8;
ceq(4) = b'*C*A*ct - 1/8;

% line 2
ceq(5) = b'*A*(c.*c) - 1/12;
ceq(6) = b'*A*(c.*ct) - 1/12;
ceq(7) = 0;%b'*A*A*c - 1/24;
ceq(8) = 0;%b'*A*A*ct - 1/24;

% line 3
ceq(9)  = b'*(ct.*c.*c) - 1/4;
ceq(10) = b'*(ct.*c.*ct) - 1/4;
ceq(11) = b'*Ct*A*c - 1/8;
ceq(12) = b'*Ct*A*ct - 1/8;

%line 4
if ~isPim2
    ceq(13) = b'*At*(c.*c) - 1/12;
    ceq(14) = b'*At*(c.*ct) - 1/12;
else
    ceq(13) = 0; %b'*At*(c.*c) - 1/12;
    ceq(14) = 0; %b'*At*(c.*ct) - 1/12;
end
ceq(15) = 0;%b'*At*A*c - 1/24;
ceq(16) = 0;%b'*At*A*ct - 1/24;

if ~isPim2
    %line 5
    ceq(17) = bt'*(c.*c.*c) - 1/4;
    ceq(18) = bt'*(c.*c.*ct) - 1/4;
    ceq(19) = bt'*C*A*c - 1/8;
    ceq(20) = bt'*C*A*ct - 1/8;
else
    ceq(17) = 0;%bt'*(c.*c.*c) - 1/4;
    ceq(18) = 0;%bt'*(c.*c.*ct) - 1/4;
    ceq(19) = 0;%bt'*C*A*c - 1/8;
    ceq(20) = 0;%bt'*C*A*ct - 1/8;
end

%line 6
ceq(21) = bt'*A*(c.*c) - 1/12;
ceq(22) = bt'*A*(c.*ct) - 1/12;
ceq(23) = 0;%bt'*A*A*c - 1/24;
ceq(24) = 0;%bt'*A*A*ct - 1/24;

%line 7
if ~isPim2
    ceq(25) = bt'*(ct.*c.*c) - 1/4;
    ceq(26) = bt'*(ct.*c.*ct) - 1/4;
    ceq(27) = bt'*Ct*A*c - 1/8;
    ceq(28) = bt'*Ct*A*ct - 1/8;
else
    ceq(25) = 0; %bt'*(ct.*c.*c) - 1/4;
    ceq(26) = 0; %bt'*(ct.*c.*ct) - 1/4;
    ceq(27) = 0; %bt'*Ct*A*c - 1/8;
    ceq(28) = 0; %bt'*Ct*A*ct - 1/8;
end

%line 8
if ~isPim2
    ceq(29) = bt'*At*(c.*c) - 1/12;
    ceq(30) = bt'*At*(c.*ct) - 1/12;
else
    ceq(29) = 0;% bt'*At*(c.*c) - 1/12;
    ceq(30) = 0;% bt'*At*(c.*ct) - 1/12;
end
ceq(31) = 0;%bt'*At*A*c - 1/24;
ceq(32) = 0;%bt'*At*A*ct - 1/24;

%line 9
ceq(33) = b'*(c.*ct.*c) - 1/4;
ceq(34) = b'*(c.*ct.*ct) - 1/4;
ceq(35) = b'*C*At*c - 1/8;
ceq(36) = b'*C*At*ct - 1/8;

%line 10
ceq(37) = b'*A*(ct.*c) - 1/12;
ceq(38) = b'*A*(ct.*ct) - 1/12;
ceq(39) = 0;%b'*A*At*c - 1/24;
ceq(40) = 0;%b'*A*At*ct - 1/24;

%line 11
ceq(41) = b'*(ct.*ct.*c) - 1/4;
ceq(42) = b'*(ct.*ct.*ct) - 1/4;
ceq(43) = b'*Ct*At*c - 1/8;
ceq(44) = b'*Ct*At*ct - 1/8;

%line 12
if ~isPim2
    ceq(45) = b'*At*(ct.*c) - 1/12;
    ceq(46) = b'*At*(ct.*ct) - 1/12;
else
    ceq(45) = 0;%b'*At*(ct.*c) - 1/12;
    ceq(46) = 0;%b'*At*(ct.*ct) - 1/12;
end
ceq(47) = 0;%b'*At*At*c - 1/24;
ceq(48) = 0;%b'*At*At*ct - 1/24;

%line 13
if ~isPim2
    ceq(49) = bt'*(c.*ct.*c) - 1/4;
    ceq(50) = bt'*(c.*ct.*ct) - 1/4;
    ceq(51) = bt'*C*At*c - 1/8;
    ceq(52) = bt'*C*At*ct - 1/8;
else
    ceq(49) = 0; %bt'*(c.*ct.*c) - 1/4;
    ceq(50) = 0; %bt'*(c.*ct.*ct) - 1/4;
    ceq(51) = 0; %bt'*C*At*c - 1/8;
    ceq(52) = 0; %bt'*C*At*ct - 1/8;
end

%line 14
ceq(53) = bt'*A*(ct.*c) - 1/12;
ceq(54) = bt'*A*(ct.*ct) - 1/12;
ceq(55) = 0;%bt'*A*At*c - 1/24;
ceq(56) = 0;%bt'*A*At*ct - 1/24;

%line 15
if ~isPim2
    ceq(57) = bt'*(ct.*ct.*c) - 1/4;
    ceq(58) = 0;% bt'*(ct.*ct.*ct) - 1/4;
    ceq(59) = bt'*Ct*At*c - 1/8;
    ceq(60) = 0;%bt'*Ct*At*ct - 1/8;
else
    ceq(57) = 0;%bt'*(ct.*ct.*c) - 1/4;
    ceq(58) = 0;% bt'*(ct.*ct.*ct) - 1/4;
    ceq(59) = 0;%bt'*Ct*At*c - 1/8;
    ceq(60) = 0;%bt'*Ct*At*ct - 1/8;
end
%line 16
if ~isPim2
    ceq(61) = bt'*At*(ct.*c) - 1/12;
else
    ceq(61) =0;% bt'*At*(ct.*c) - 1/12;
end
ceq(62) = 0;%bt'*At*(ct.*ct) - 1/12;
ceq(63) = 0;%bt'*At*At*c - 1/24;
ceq(64) = 0;%bt'*At*At*ct - 1/24;

end
