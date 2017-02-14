function [X] = pack_implicit(A,b)
% File: pack_implicit.m
% Author: Sidafa Conde
% Email: sconde@umassd.edu
% School: UMass Dartmouth
% Date: 02/01/2017
% Purpose: packs (column major) the Butcher coefficient

s = length(b);
assert(isequal(size(A,1),s) && isequal(size(A,2),s),...
'dimension mismatched');

ind = tril(true(s));
x_a = A(ind);
X = [x_a; b(:)];
end
