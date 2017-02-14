function [ceq] = nlc_non_imex(x, s, k,pex, pim, plin, implicit_type, special_assumption, enforce_positivity)

if strcmpi(special_assumption, 'g')
    [con, ceq] = nlc_imex_g(x, s, k,pex, pim, plin, implicit_type, special_assumption, enforce_positivity);
elseif strcmpi(special_assumption, 'c')
    [con, ceq] = nlc_imex_c(x, s, k,pex, pim, plin, implicit_type, special_assumption, enforce_positivity);
else
    error('nlc-special-assumption','type not recognized');
end
con;

end
