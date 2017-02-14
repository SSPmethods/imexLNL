function [A, b, Ahat, bhat] = shuOsher2Butcher(alpha, alpha_hat, c, ct, s)
%function [Re,P,Q] = shuOsher2Butcher(alpha, alpha_hat, c, ct, s)
%Converting Modified Shu-Osher form to Butcher

LAMBDA = alpha;
LAMBDA_HAT = alpha_hat;

I = speye(size(LAMBDA));

K = (I - (LAMBDA + LAMBDA_HAT))\(LAMBDA/c);
KT = (I - (LAMBDA + LAMBDA_HAT))\(LAMBDA_HAT/ct);

A = K(1:s,1:s); b = K(s+1,1:s);
Ahat = KT(1:s, 1:s); bhat = KT(s+1,1:s);

end