function [A, b, c, At, bt, ct, r1, r2] = unpack_imex(x, s, k, implicit_type, special_assumption)

% SDIRK-G               : x = [A(:), b(:), At(:), bt(:), gam, r1]
% {SDIRK-B, SDIRK-BC}   : x = [A(:), b(:), At(:), gam, r1]
% DIRK-G                : x = [A(:), b(:), At(:), bt(:), r1]
% {DIRK-B, DIRK-BC}     : x = [A(:), b(:), At(:), r1]
% DIRK-BC-II            : x = [A(:), b(:), At(1:s-1:), r1]

xx = x(1:end-1);
if ~isinf(k)
    % if k is not infinity, then optimize over r (the explicit ssp coef)
    r1 = -x(end);
    r2 = k*r1;
else
    % for k = inf, then we revers the optimization variable to implcit
    % this way, r (explicit) = 0
    r2 = -x(end);
    r1 = 0;
end

if strcmpi(implicit_type, 'sdirk')
    gam = xx(end);
    xx = xx(1:end-1);
end

n_exp = 0.5*s*(s-1) + s;
xx_exp = xx(1:n_exp);
xx_imp = xx(n_exp+1:end);

% explicit method
A = zeros(s);
xx_a = xx_exp(1:end-s);
A(tril(true(s),-1)) = xx_a;
b = xx_exp(end-s+1:end);

% implicit method
At = zeros(s);

if strcmpi(implicit_type, 'dirk')
    if strcmpi(special_assumption, 'c')
        at_n = 0.5*s*(s-1) + (s-1) -2;
        ind = tril(true(s)); ind(1,1) = false; ind(2,:) = false;
        xx_at = xx_imp(1:at_n); a21_half = A(2,1)/2;
        At(ind) = xx_at; At(2,[1 2]) = a21_half;
        
    elseif strcmpi(special_assumption,'G')
        at_n = 0.5*s*(s+1);
        xx_at = xx_imp(1:at_n);
        At(tril(true(s))) = xx_at;
    else
        error('unpack-special-assumption','type not recognized');
    end
elseif strcmpi(implicit_type, 'sdirk')
    at_n = 0.5*s*(s-1);
    xx_at = xx_imp(1:at_n);
    At(tril(true(s),-1)) = xx_at;
    At = At + gam*eye(size(At));
end

if strcmpi(special_assumption,'b')
    bt = b;
else
    bt = xx_imp(at_n +1: at_n+s);
end


% k = [xx_a; xx_at];
% %warning off;
% options = optimoptions('fsolve','Display','off','TolFun',1e-16, 'TolX',1e-16);
% [K, fval] = fsolve(@(k) enforce_c(k,s),k, options);
% a_n = 0.5*s*(s-1);
% at_n = a_n + (s-1);
%
% x_a = K(1:a_n);
% x_at = K(a_n + 1:end); assert(isequal(numel(x_at), at_n));
%
% A = zeros(s); A(tril(true(s),-1)) = a_n;
%
% at_n = 0.5*s*(s-1) + (s-1);
% ind = tril(true(s)); ind(1,1) = false;
% At(ind) = x_at;

c = sum(A,2); ct = sum(At, 2);
end


