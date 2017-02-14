function [con, ceq] = nlc_imex(x, s, k,pex, pim, plin, implicit_type, special_assumption, enforce_positivity, useSameBSameC, varargin)

if useSameBSameC
    % only the fourth order method at this time
%     assert(isequal(pex, 4) && isequal(pim, 4) && isequal(plin, 4),...
%         'only fourth order method at this time');
    
    [con, ceq] = nlc_imex_sameB_sameC(x, s, k,pex, pim, plin, implicit_type, special_assumption, enforce_positivity);

else
    
    if isempty(varargin) % regular order condition
        % making no distinction between the explicit and implicit linear order
        % conditions (i.e P^e_{lin} = P^i_{lin} = plin
        [con, ceq] = nlc_imex_g(x, s, k,pex, pim, plin, implicit_type, special_assumption, enforce_positivity);
    else
        % making distinction between the linear order conditions
        pi_lin = varargin{1};
        [con, ceq] = nlc_imex_g_diff_plin(x, s, k,pex, pim, plin, implicit_type, special_assumption, enforce_positivity, pi_lin);
    end
end
end
