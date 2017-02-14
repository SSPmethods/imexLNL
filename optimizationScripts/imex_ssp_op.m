% restart = 0;  random start
% restart = 1;  start from known method
% restart = 2;  start with perturbation of loaded method
clear all; close all; clc
addpath('utils/'); addpath('order_cond/'); addpath('helpers/');

nonSSP_search = false;
am_minref = rand;
implicit_type = 'DIRK';
special_assumption = 'G';
enforce_positivity = 'all';
displayType = 'iter'; % off/iter
parallelOption = 'never'; % Never/ Always
useSameBSameC = false;

s = 3;
pex = 3;
pim = 3;
plin = 3;
k = 0;
restart = 0;
minr = 0.1;
minr_bound = s;

starttype = 'random';
rng_state = rng('shuffle');
[n_exp, n_imp] = set_n_imex(s,implicit_type, special_assumption);
n = n_exp + n_imp;

if restart ==0
    x0 = rand(n,1); x0(end) = 0;
    r = 0;
elseif restart == 1
    x0 = X;
elseif restart == 2
    x0(1:end-1) = X(1:end-1) + rand(size(X(1:end-1))); x0(end) = -minr - 0.1;
end

% The conjectured ssp bound for explicit methods is at C <= s
lb = -5 + zeros(1,n); lb(end) = -minr_bound;
ub = 5 + zeros(1,n); ub(end) = -minr;

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

nlc_imex_func = @nlc_imex;

warning('off');

while ((r < minr) || (info == -2) || info == 0) && (count < 10) && ~(ACCURACY <= -15)

    [X, FVAL, info, output] = fmincon(@am_obj, x0, [],[], Aeq, beq, lb,...
        ub,@(x) nlc_imex_func(x, s, k,pex, pim, plin, implicit_type, special_assumption, enforce_positivity, false), opts);

    r = -FVAL; if isempty(r); r = abs(X(end));end
    count = count + 1;

    % check the accuracy of the method
    [~, ceq] = nlc_imex_func(X, s, k,pex, pim, plin, implicit_type, special_assumption, enforce_positivity, false);
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

[A, b, c, At, bt, ct, r, rt] = unpack_imex(X, s, k, implicit_type, special_assumption);
[v, alpha, alpha_hat]  = Butcher2ShuOsher(A, At, b', bt', r, rt);
[con, ceq] = nlc_imex(X, s, k,pex, pim, plin, implicit_type, special_assumption, enforce_positivity, false);
ACCURACY = -floor(min(real(-log10(ceq))));
saveCondition = (ACCURACY <= -14);

positiveCoeff = ~any(con > 1e-13);
saveCondition = saveCondition && positiveCoeff;

fprintf(1, '\t%s Accuracy = %d, positiveCoeff = %d\n', repmat('*',1,7),...
    ACCURACY, positiveCoeff);

