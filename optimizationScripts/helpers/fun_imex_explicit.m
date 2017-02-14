function [ ceq] = fun_imex_explicit(x, A, b, r, k,pex, pim, plin, implicit_type, special_assumption, enforce_positivity, useSameBSameC, varargin)
% File: fun_imex_explicit.m
% Author: Sidafa Conde
% Email: sconde@umassd.edu
% School: UMass Dartmouth
% Date: 02/01/2017
% Purpose: %TODO

[At, bt, ~] = unpack_rk(x, numel(b), implicit_type);

c = sum(A,2); ct = sum(At, 2);
b = b(:); bt = bt(:); Ct = diag(ct);

con = ssp_cond_imex(numel(b), A, b, r, At, bt, k, false);

if (pex >=2 && pim >= 2)
    % First order condition
    ceq(1) = sum(bt) - 1.0;
    
    % Second order conditons (all linear orders)
    ceq(2) = b'*ct - 1/2;
    ceq(3) = bt'*c - 1/2;
    ceq(4) = bt'*ct - 1/2;
end

if (pex >=3 )
    % explicit Nonlinear 3rd order    
    ceqtemp = p3_nonlinear_couple( A, b, c, At, bt, ct);
    ceq = [ceq(:); ceqtemp(:)];
end

if (pim >= 3)
    % implicit Nonlinear 3rd order
    ceq(end+1) = bt'*(ct.^2) - 1/3;
end

if (pex == 4 )
    % explicit Nonlinear 4th order
    ceq_temp4 = p4_nonlinear_couple( A, b, c, At, bt, ct);
    ceq = [ceq(:); ceq_temp4(:)];
end

if (pim >= 4)
    % implicit Nonlinear 4th order
    ceq(end+1) = bt'*(ct.^3) - 1/4;
    ceq(end+1) = bt'*Ct*At*ct - 1/8;
    ceq(end+1) = bt'*At*(ct.*ct) - 1/12;
end

 %enforce no-zero in diagonal
%diagAt = min(abs(diag(At)));
%diagAt = min(abs(diag(At)));
%diagAtCond = diagAt > 1e-12;

%diagAtCond = 0.1 - diagAtCond;
%con = [con(:); diagAtCond];


% Coupled-linear-order-condition
ceq_lin = general_coupled_linear(plin, A, b, At, bt);
ceq  = [ceq(:); ceq_lin(:)];

if strcmpi(implicit_type, 'sdirk')
    % extract gamma
    gam = diag(At);
    assert(~any(diff(gam)),'gamma-check', 'method is not really SDIRK');
    gam = gam(2); % pick the second one
    con = [con(:); -gam];
end
end
