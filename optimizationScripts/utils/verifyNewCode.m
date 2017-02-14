%verifyNewCode.m
clear all; close all;clc

load r1_2.500000_r2_2.200000e+00_am_r_5.500000.mat
type = 'B';

[A_n,b_n,c_n, At_n, bt_n, ct_n, r1_n, r2_n] = unpack_imex(X, s, type);

%accuracy of correct unpack
[con, ceq] = nlc_imex_typeB(s, p, A, b, c, At, bt, ct, r1, r2);

%accuracy of new unpack
% this unpack is not the same as the old unpack
[con_n, ceq_n] = nlc_imex_typeB(s, p, A_n, b_n, c_n, At_n, bt_n, ct_n, r1_n, r2_n);

assert(norm(A - A_n) < 1e-16,'As are not the same');
assert(norm(b - b_n) < 1e-16,'As are not the same');
assert(norm(c - c_n) < 1e-16,'As are not the same');