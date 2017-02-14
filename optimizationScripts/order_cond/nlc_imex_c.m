function [con, ceq] = nlc_imex_c(x, s, k,pex, pim, plin, implicit_type, special_assumption)
% function [con, ceq] = nlc_imex_c(x, s, k,pex, pim, plin, implicit_type, special_assumption)
% Purpose: order condition for ARK IMEX method
% 		with the assumption that b = bt and c = ct

[A, b, c, At, bt, ct, r1, r2] = unpack_imex(x, s, k, implicit_type, special_assumption);

con = ssp_cond_imex(s, A, b, r1, At, bt, r2);
b = b(:); bt = bt(:);


if (pex == 4 && pim == 3 && plin == 4)
    [ceq] = pex4_pex_lin4_pim3_pim_lin4(plin, A, b, c, At, bt, ct);
elseif (pex == 4 && pim == 3 && plin >= 5)
    % gets all the higher order linear coupled linear condition
    % starting at plin = 5
    [ceq] = pex4_pex_lin5_pim3_pim_lin5(plin, A, b, c, At, bt, ct);
else
    if (pex >=2 && pim >= 2)
        % First order condition
        ceq(1) = sum(b) - 1.0;
        ceq(2) = sum(bt) - 1.0;
        
        % Second order conditons -- all coupled order condition
        ceq(3) = b'*c - 1/2;
        ceq(4) = b'*ct - 1/2;
        ceq(5) = bt'*c - 1/2;
        ceq(6) = bt'*ct - 1/2;
    end
    
    if (pex >=3 && pim >= 3)
        % Third order conditions
        ceq(end+1) = b'*c.^2 - 1/3;
        ceq(end+1) = bt'*(ct.^2) - 1/3;
        ceqtemp = p3_nonlinear_couple( A, b, c, At, bt, ct);
        ceq = [ceq(:); ceqtemp(:)];
    end
    
    if (pex == 4 )
        ceq_temp4 = p4_nonlinear_couple( A, b, c, At, bt, ct);
        ceq = [ceq(:); ceq_temp4(:)];
    end
    
    % Coupled-linear-order-condition
    ceq_lin = general_coupled_linear(plin, A, b, At, bt);
    ceq  = [ceq(:); ceq_lin(:)];
end
ceq_b = b - bt; % same b : equal weight
ceq_c = c - ct; % same c : equal abscissas (canopy)
ceq = [ceq(:); ceq_b(:); ceq_c(:)];


end
