function [con, ceq] = nlc_so_imex_g(x, s, k,pex, pim, plin)
% function [con, ceq] = nlc_so_imex_g(x, s, k,pex, pim, plin)
% Shu-Osher optimization

[A, b, c, At, bt, ct, r1, r2] = unpack_so_imex(x, s, k);

con = ssp_cond_imex(s, A, b, r1, At, bt, r2);
    % the first two-order conditions
    b = b(:); bt = bt(:); %c = sum(A,2); ct = sum(At, 2);
if (pex == 4 && pim == 3 && plin == 4)
    [ceq] = pex4_pex_lin4_pim3_pim_lin4(plin, A, b, c, At, bt, ct);
elseif (pex == 4 && pim == 3 && plin == 5)
    [ceq] = pex4_pex_lin5_pim3_pim_lin5(plin, A, b, c, At, bt, ct);
elseif (pex == 4 && pim == 3 && plin == 6)
    [ceq] = pex4_pex_lin6_pim3_pim_lin6(plin, A, b, c, At, bt, ct);
else
    %first order condition (all verified)
    ceq(1) = sum(b) - 1.0;
    ceq(2) = sum(bt) - 1.0;
    
    %second order conditons ( all verified)
    ceq(3) = b'*c - 1/2;
    ceq(4) = b'*ct - 1/2;
    ceq(5) = bt'*c - 1/2;
    ceq(6) = bt'*ct - 1/2;
    
    if pex >= 3
        ceq(end+1) = b'*c.^2 - 1/3;
        ceqtemp = p3_nonlinear_couple( A, b, c, At, bt, ct);
        ceq = [ceq(:); ceqtemp(:)];
    end
    
    if pim >= 3
        ceq(end+1) = bt'*(ct.^2) - 1/3;
    end
    
    if pex >= 4
       ceq_temp4 = p4_nonlinear_couple( A, b, c, At, bt, ct);
       ceq = [ceq(:); ceq_temp4(:)];
    end

    ceq_lin = general_coupled_linear(plin, A, b, At, bt);
    ceq  = [ceq(:); ceq_lin(:)];
 end
end
