function [n_exp, n_imp] = set_n_imex(s,implicit_type, special_assumption, fsolve_c_ct)

a_n = 0.5*s*(s-1);
b_n = s;
r1 = 1;

if strcmpi(implicit_type, 'dirk')
    % implicit-type
    switch lower(special_assumption)
        case 'g'
            at_n = s;
        case 'c' % case with abscissa being equal
            at_n = s - 1;
        otherwise
            at_n = 0.5*s*(s+1) -2;
    end
elseif strcmpi(implicit_type, 'sdirk')
    at_n = 1; %just gamma
end

at_n = at_n + a_n;
% assumption
switch lower(special_assumption)
    case {'bc','b'}
        bt_n = 0;
    otherwise
        bt_n = s;
end

n_exp = a_n + b_n + r1;
n_imp = at_n + bt_n;

end
