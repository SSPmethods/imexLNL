s = numel(b);
[A, b, c, At, bt, ct, r, rt] = unpack_imex(X, s, k, implicit_type, special_assumption);

[con, ceq] = nlc_imex(X, s, k,pex, pim, plin, implicit_type, special_assumption, enforce_positivity);
ACCURACY = -floor(min(real(-log10(ceq))));
saveCondition = (ACCURACY <= -14);

saveCondition = saveCondition && ~any(con > 1e-15);
repeat = 0;
restart = 1;
%alwasy assumer first that I can start pushing
saveCondition = 1;
minr_old = minr;
pusher = 1;
enforce_positivity = 1; free_optimize =1;

while (saveCondition) && (repeat < 5)
    
    duplicate_method = 0;
    XX = X;
    imex_ssp_op;
    
    
    [A, b, c, At, bt, ct, r, rt] = unpack_imex(X, s, k, implicit_type, special_assumption);
    [con, ceq] = nlc_imex(X, s, k,pex, pim, plin, implicit_type, special_assumption, enforce_positivity);
    ACCURACY = -floor(min(real(-log10(ceq))));
    
    
    if repeat >= 10
        break
    end
    
    if (~saveCondition) %|| (duplicate_method)
        X = XX;
        saveCondition = 1;
        repeat = repeat + 1
        pusher = 0.5*pusher;
    else
        repeat = 0
        pusher = 1;
        restart = 1;
        minr_old = r;
    end
    
    minr = minr_old + pusher;
    minr = r;
end
