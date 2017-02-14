% match_zackBestMethod.m
% Purpose: try to see if the optimizer will converge to the optimal method
% Zack got

clear all; close all; clc

addpath('utils/');
addpath('order_cond/');
addpath('helpers/');
addpath('optSrc/');

zack = 0;
displayType = 'iter';
restart = 1;
enforce_positivity = 'all';
free_optimize =0;

if zack
    zk_rk = load('/Users/sconde/Downloads/SSPIMEXRKs_7pex_4pim_4plin_6K_.mat');
    
    
    X = AB_to_X(zk_rk.A, zk_rk.b, zk_rk.Ahat, zk_rk.bhat, zk_rk.r);
    
    
    implicit_type = 'DIRK';
    special_assumption = 'G';
    s = zk_rk.s;
    r = zk_rk.r;
    pex = zk_rk.pex;
    pim = zk_rk.pim;
    plin = zk_rk.plin;
    k = 0;
    %minr = zk_rk.r;
    %minr_bound = zk_rk.r + 1e-4;
    
else
    load('Method/DIRK/G/Pex4/Pim4/Plin6/S7/K0/method_typeG_r1_1.5000000000000_acc_-15.mat');
end

minr = 1.6;%1.826937531629875;
minr_bound = 1.6;%1.826937531629875;
n = length(X);
k = k;
imex_ssp_op
