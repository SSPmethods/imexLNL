function isNearEqual = isequalVectorEps(b, bt)
% function isqual = isequalVectorEps(b, bt)
% Purpose: determine is two vectors are equal to floating point precision

% first clean up
b(abs(b) < eps) = 0;
bt(abs(bt) < eps) = 0;

isNearEqual = all(abs(b-bt) < 1e4*eps(min(abs(b),abs(bt))));

end