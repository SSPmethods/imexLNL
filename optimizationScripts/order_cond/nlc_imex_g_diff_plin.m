function [con, ceq] = nlc_imex_g_diff_plin(x, s, k,pex, pim, plin, implicit_type, special_assumption, enforce_positivity, pi_lin)

[A, b, c, At, bt, ct, r1, r2] = unpack_imex(x, s, k, implicit_type, special_assumption);

con = ssp_cond_imex(s, A, b, r1, At, bt, r2, enforce_positivity);
b = b(:); bt = bt(:); C = diag(c); Ct = diag(ct);


if (pex >=2 && pim >= 2)
    % First order condition
    ceq(1) = sum(b) - 1.0;
    ceq(2) = sum(bt) - 1.0;
    
    % Second order conditons (all linear orders)
    ceq(3) = b'*c - 1/2;
    ceq(4) = b'*ct - 1/2;
    ceq(5) = bt'*c - 1/2;
    ceq(6) = bt'*ct - 1/2;
end



% for now, just looking at second order nonlinear order contions
% these will be just the tall-tree conditions


% phat_lin = (pi_lin - 2);
% 
e = ones(s,1);
% ceqfun = @(temp) arrayfun(@(i) temp*At^(i-1)*e - 1/factorial(i+2), 1:phat_lin);
% 
% keyboard
% ceq(7) = ceqfun(b'*A);
% ceq(8) = ceqfun(b'*At);
% ceq(9) = ceqfun(bt'*A);
% ceq(10) = ceqfun(bt'*At);

if pi_lin >= 3
    ceq(7)  = b'*A*At*e    - 1/6;
    ceq(8)  = b'*At*At*e   - 1/6;
    ceq(9)  = bt'*A*At*e   - 1/6;
    ceq(10) = bt'*At*At*e  - 1/6;
end

if pi_lin >= 4
    ceq(7)  = b'*A*At*At*e    - 1/factorial(pi_lin);
    ceq(8)  = b'*At*At*At*e   - 1/factorial(pi_lin);
    ceq(9)  = bt'*A*At*At*e   - 1/factorial(pi_lin);
    ceq(10) = bt'*At*At*At*e  - 1/factorial(pi_lin);
end

if pi_lin >= 5
    ceq(7)  = b'*A*At*At*At*e    - 1/factorial(pi_lin);
    ceq(8)  = b'*At*At*At*At*e   - 1/factorial(pi_lin);
    ceq(9)  = bt'*A*At*At*At*e   - 1/factorial(pi_lin);
    ceq(10) = bt'*At*At*At*At*e  - 1/factorial(pi_lin);
end

end
