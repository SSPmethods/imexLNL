function [A, b] = unpack_implicit(X,s)
% File: unpack_implicit.m
% Author: Sidafa Conde
% Email: sconde@umassd.edu
% School: UMass Dartmouth
% Date: 02/01/2017
% Purpose: unpacks (column major) into Butcher coefficient

n_at = 0.5*s*(s+1);
ind = tril(true(s));
A = zeros(s);
A(ind) = X(1:n_at);
b = X(n_at+1:end);
end
