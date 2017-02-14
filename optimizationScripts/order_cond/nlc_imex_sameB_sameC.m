function [con, ceq] = nlc_imex_sameB_sameC(x, s, k,pex, pim, plin, implicit_type, special_assumption, enforce_positivity)

[A, b, c, At, bt, ct, r1, r2] = unpack_imex(x, s, k, implicit_type, special_assumption);

con = ssp_cond_imex(s, A, b, r1, At, bt, r2, enforce_positivity);
b = b(:); bt = bt(:); C = diag(c); Ct = diag(ct);


assert(isequal(pex, pim),...
    'Only pex == pim right now');
b = b(:); bt = bt(:);

%con = ssp_cond_imex(s, A, b, ssp104r, At, bt, rt, true);

c = sum(A,2); ct = sum(At,2);

p = pex;

% conditions for b = bt and c = ct;
ceq_sp = b - bt;
ceq_sp = [ceq_sp; (c - ct)];

if p >= 1
    ceq(1) = sum(b) - 1;
end

if p >= 2
    ceq(2) = b'*c - 1/2;
end

if p >= 3
    ceq(3) = b'*A*c - 1/6;
    ceq(4) = b'*At*c -1/6;
    ceq(5) = b'*c.^2 - 1/3;
end

if p == 4
    ceq(6) = b'*c.^3 - 1/4;
    ceq(7) = b'*C*A*c - 1/8;
    ceq(8) = b'*C*At*c - 1/8;
    ceq(9) = b'*A*c.^2 - 1/12;
    ceq(10) = b'*A*A*c - 1/24;
    ceq(11) = b'*A*At*c - 1/24;
    ceq(12) = b'*At*c.^2 - 1/12;
    ceq(13) = b'*At*A*c - 1/24;
    ceq(14) = b'*At*At*c - 1/24;
end

ceq = [ceq(:); ceq_sp(:)];



% Coupled-linear-order-condition
ceq_lin = general_coupled_linear(plin, A, b, At, bt);
ceq  = [ceq(:); ceq_lin(:)];
end
