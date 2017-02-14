function X = AB_to_X(A, b, Ahat, bhat, r)
% Sidafa Conde & Sigal Gottlieb
% SSP IMEX
% 03/23/2015
% given A, b, Ahat, bhat, r
% pack to X (optimization variable)
% infers stage number (s) by length of b and bhat

assert(isequal(length(b), length(bhat)),'b-bhat-length',...
    '(inferring) stage number from b and bhat are unequal');

s = length(b);
b = b(:); bhat = bhat(:);
exp_ind = tril(true(s),-1);
imp_ind = tril(true(s));
xx_exp = A(exp_ind);
xx_imp = Ahat(imp_ind);
X = [xx_exp; b; xx_imp; bhat; -abs(r)];


end
