% Purpose: given an optimal explicit method, solve the order condition for
% the corresponding implicit method ( this coincides with k = 0 )

addpath('utils/'); addpath('order_cond/');

%restart = 0;
implicit_type = 'DIRK';
special_assumption = 'G';
tol = 1e-15;

k = 0;
starttype = 'random';
if restart ==0
    [~, n_imp] = set_n_imex(s,implicit_type, special_assumption);
    n = n_imp;
    rng_state = rng('shuffle');
    x0 = rand(( n_imp),1); %x0(end) = 0;
elseif restart == 1
    x0 = X; n = length(x0);
elseif restart == 2
    x0(1:end-1) = X(1:end-1) + rand(size(X(1:end-1))); x0(end) = -minr - 0.1;
end


% first try to solve the nonlinear order conditions
fsolve_option = optimoptions('fsolve','Display','iter','TolFun',1e-14, 'MaxIter',500000);
fun = @(x) fun_imex_explicit(x, A, b, r, k, pex, pim, plin, implicit_type, special_assumption, enforce_positivity, useSameBSameC);
[x0, fval] = fsolve(fun, x0);

[At, bt] = unpack_implicit(x0,s);

keyboard



% The conjectured ssp bound for explicit methods is at C <= s
lb = -1 + zeros(1,n); %lb(end) = -minr_bound;
ub = 1 + zeros(1,n); %ub(end) = -minr;

%TODO: should find the index for the implicit method
% and make sure they are strickly positive

tol = 1e-15;
opts=optimset('MaxFunEvals',8500,'TolCon',tol,'TolFun',tol,'TolX',tol,...
    'GradObj','on','MaxIter',2000,'Diagnostics','off','Display',displayType,...
    'UseParallel',parallelOption,'Algorithm','sqp');

Aeq = [];
beq = [];

info = -2; FVAL = 0;
am_r = 0; count = 0;
forceStart = 0;
ACCURACY = 5;

if useSameBSameC
    nlc_imex_func = @nlc_imex_reduced;
else
    nlc_imex_func = @nlc_imex_explicit;
end

displayType = 'iter';
opts = optimset('MaxFunEvals',10000,'TolCon',tol,'TolFun',tol,'TolX',tol,...
    'GradObj','on','MaxIter',5000,'Diagnostics','off','Display',displayType,...
    'UseParallel',parallelOption,'Algorithm','sqp');

saveCondition = 0;
count = 0;
Aeq = [];
beq = [];

info = -2; FVAL = 0;
am_r = 0; count = 0;
forceStart = 0;
ACCURACY = 5;

while ((info == -2) || info == 0) && (count < 10) && ~(ACCURACY <= -15)
    
    [X, FVAL, info, output] = fmincon(@(x) am_obj(x, r), x0, [],[], Aeq, beq, lb,...
        ub,@(x) nlc_imex_explicit(x, A, b, r, k, pex, pim, plin, implicit_type, special_assumption, enforce_positivity, useSameBSameC), opts);
    
    r = -FVAL;
    if isempty(r)
        r = abs(X(end));
    end
    
    count = count + 1;
    
        
    % check the accuracy of the method
    [~, ceq] = nlc_imex_explicit(X, A, b, r, k,pex, pim, plin, implicit_type, special_assumption, enforce_positivity, useSameBSameC);
    ACCURACY = -floor(min(real(-log10(ceq))));
    
    
    if ACCURACY <= -10
        x0 = X;
    else
        x0(1:end-1) = X(1:end-1) + 0.1*rand(size(X(1:end-1)));
    end
    
    if floor(log10(output.constrviolation)) <= -15
        break
    end
    
end

[At, bt, ~] = unpack_rk(X, numel(b), implicit_type);
X = AB_to_X(A, b, At, bt, r);

[A, b, c, At, bt, ct, r, rt] = unpack_imex(X, s, k, implicit_type, special_assumption);
[v, alpha, alpha_hat]  = Butcher2ShuOsher(A, At, b', bt', r, rt);



[con, ceq] = nlc_imex(X, s, k,pex, pim, plin, implicit_type, special_assumption, false, false);
ACCURACY = -floor(min(real(-log10(ceq))));
saveCondition = (ACCURACY <= -14)

% RKStabiltyPlotLeveque(At,bt(:)',-100,100,-100,100,{'b'}, true);
% keyboard

positiveCoeff = ~any(con > 1e-13);
saveCondition = saveCondition && positiveCoeff;

fprintf(1, '\t%s Accuracy = %d, positiveCoeff = %d\n', repmat('*',1,7),...
    ACCURACY, positiveCoeff);

    keyboard

check_and_save_method;
