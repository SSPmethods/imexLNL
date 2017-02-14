%clear all; close all; clc
% 
% fil_method = dir('~/Desktop/*.mat');
% 
% for i = 1:numel(fil_method)
%     load(['~/Desktop/' fil_method(i).name]);
    
    implicit_type = 'DIRK';
    special_assumption ='G';    % G/B/BC
    %r = abs(r1); rt = abs(r2);
    %k = rt/r;
    % explicit method
    A = [A;0.1*zeros(1,s)];
    A = [A zeros(s+1,1)];
    
    At = [At;0.1*zeros(1,s)];
    At = [At zeros(s+1,1)];
    b = [b; 0];
    bt = [bt; 0];
    s = size(A,1);
    
    exp_ind = tril(true(s),-1);
    xx_exp = A(exp_ind);
    
    %implicit method
    imp_ind = tril(true(s), 0);
    xx_imp = At(imp_ind);
    
    X = [xx_exp; b(:); xx_imp; bt(:); -abs(r)];
    
    [A, b, c, At, bt, ct, r, rt] = unpack_imex(X, numel(b), k, implicit_type, special_assumption);
   
    
    [con, ceq] = nlc_imex(X, s, k,  pex, pex_lin, pim, pim_lin, implicit_type, special_assumption);
    ACCURACY = -floor(min(real(-log10(ceq))))
    saveCondition = (ACCURACY <= -14);
    
    saveCondition = saveCondition && ~any(con > 1e-15)
    
    
    if saveCondition
        [A, b, c, At, bt, ct, r, rt] = unpack_imex(X, s, k, implicit_type, special_assumption);
        [v, alpha, alpha_hat]  =Butcher2ShuOsher(A, At, b', bt', r, k);
        [con, ceq] = nlc_imex(X, s, k,  pex, pex_lin, pim, pim_lin, implicit_type, special_assumption);
        ACCURACY = -floor(min(real(-log10(ceq))));
        
        check_and_save_method;
    end
