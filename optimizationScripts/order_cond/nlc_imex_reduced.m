function [con, ceq] = nlc_imex_reduced(x, s, k,pex, pim, plin, implicit_type, special_assumption, enforce_positivity, useSameBSameC)

% first unpack the IMEX-RK method
[A, b, c, At, bt, ct, r, rt] = unpack_imex(x, s, k, implicit_type, special_assumption);

b = b(:); bt = bt(:);

con = ssp_cond_imex(s, A, b, r, At, bt, rt, true);

C = diag(c);
Ct = diag(ct);

% conditions for b = bt and c = ct;
ceq_sp = b - bt;
ceq_sp = [ceq_sp; (c - ct)];

ceq(1) = sum(b) - 1;

% pex = pim = 2;
if (pex >= 2) && (pim >= 2)
    ceq(2) = b'*c - 1/2;
end

if (pex >= 3) 
    % pex = pim = 3
    ceq(end+1) = b'*c.^2 - 1/3;
    
    % plin = 3
    ceq(end+1) = b'*A*c - 1/6;
    ceq(end+1) = b'*At*c -1/6;
end

if (pex == 4)
    ceq(end+1) = b'*c.^3 - 1/4;
    ceq(end+1) = b'*C*A*c - 1/8;
    ceq(end+1) = b'*A*c.^2 - 1/12;
end

if (pim == 4)
    ceq(end+1) = b'*C*At*c - 1/8;
    ceq(end+1) = b'*At*c.^2 - 1/12;
end

if (plin >= 4)
    ceq(end+1) = b'*A*A*c - 1/24;
    ceq(end+1) = b'*A*At*c - 1/24;
    ceq(end+1) = b'*At*A*c - 1/24;
    ceq(end+1) = b'*At*At*c - 1/24;
end

if (plin >= 5)
   ceq(end+1) = b'*A*A*A*c - 1/120;
   ceq(end+1) = b'*At*A*A*c - 1/120;
   ceq(end+1) = b'*At*At*A*c - 1/120;
   ceq(end+1) = b'*At*At*At*c - 1/120;
   ceq(end+1) = b'*A*At*A*c - 1/120;
   ceq(end+1) = b'*A*At*At*c - 1/120;
   ceq(end+1) = b'*A*A*At*c - 1/120;
end

if (plin >= 6)
    ceq(end+1) = b'*A*A*A*A*c - 1/720;
    ceq(end+1) = b'*At*A*A*A*c - 1/720;
    ceq(end+1) = b'*At*At*A*A*c - 1/720;
    ceq(end+1) = b'*At*At*At*A*c - 1/720;
    ceq(end+1) = b'*At*At*At*At*c - 1/720;
    ceq(end+1) = b'*A*At*A*A*c - 1/720;
    ceq(end+1) = b'*A*At*At*A*c - 1/720;
    ceq(end+1) = b'*A*At*At*At*c - 1/720;
    ceq(end+1) = b'*A*A*At*A*c - 1/720;
    ceq(end+1) = b'*A*A*At*At*c - 1/720;
    ceq(end+1) = b'*A*A*A*At*c - 1/720;
end

ceq = [ceq(:); ceq_sp(:)];

end
