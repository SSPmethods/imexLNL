function [con, ceq] = nlc_imex_g(x, s, k,pex, pim, plin, implicit_type, special_assumption, enforce_positivity)
% File: nlc_imex_g.m
% Author: Sidafa Conde
% Email: sconde@umassd.edu
% School: UMass Dartmouth
% Last modified: 02/07/2017
% Purpose:  %TODO

[A, b, c, At, bt, ct, r1, r2] = unpack_imex(x, s, k, implicit_type, special_assumption);

con = ssp_cond_imex(s, A, b, r1, At, bt, r2, enforce_positivity);
b = b(:); bt = bt(:); C = diag(c); Ct = diag(ct);


if (pex >=2 && pim >= 2)
    % First order condition
    ceq(1) = sum(b) - 1.0;
    ceq(2) = sum(bt) - 1.0;

    % Second order conditons (all linear orders, i.e. tall tree)
    ceq(3) = b'*c - 1/2;
    ceq(4) = b'*ct - 1/2;
    ceq(5) = bt'*c - 1/2;
    ceq(6) = bt'*ct - 1/2;
end

if (pex >=3 )
    % explicit Nonlinear 3rd order
    ceq(end+1) = b'*c.^2 - 1/3;
    ceqtemp = p3_nonlinear_couple( A, b, c, At, bt, ct, isequal(pim, 2));
    ceq = [ceq(:); ceqtemp(:)];
    
end

if (pim >= 3)
    % implicit Nonlinear 3rd order
    ceq(end+1) = bt'*(ct.^2) - 1/3;
end

if (pex == 4 )
    % explicit Nonlinear 4th order
    ceq(end+1) = b'*c.^3 - 1/4;
    ceq(end+1) = b'*C*A*c - 1/8;
    ceq(end+1) = b'*A*(c.*c) - 1/12;

    % TODO: need to differentiate the pim = 3 case
    ceq_temp4 = p4_nonlinear_couple( A, b, c, At, bt, ct, isequal(pim, 2));
    ceq = [ceq(:); ceq_temp4(:)];
end

if (pim >= 4)
    % implicit Nonlinear 4th order

    ceq(end+1) = bt'*(ct.^3) - 1/4;
    ceq(end+1) = bt'*Ct*At*ct - 1/8;
    ceq(end+1) = bt'*At*(ct.*ct) - 1/12;
end

% Coupled-linear-order-condition
ceq_lin = general_coupled_linear(plin, A, b, At, bt);
ceq  = [ceq(:); ceq_lin(:)];

end
