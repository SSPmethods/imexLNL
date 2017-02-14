clear all; %close all; clc
addpath('order_cond/'); 
format
pex = 2;
pim = 2;
k = 0;
implicit_type = 'DIRK';
special_assumption = 'G';
enforce_positivity = 0;

S = 2:10;
P = 2:10;

C = zeros(numel(S)+1,numel(P)+1);
C(1,2:end) = P; C(2:end,1) = S;

for j = 1:numel(P)
    p = P(j);
    for i = 1:numel(S)
        s = S(i);
        try
            direc = sprintf('Method/%s/%s/Pex%d/Pim%d/Plin%d/S%d/K%s/',...
                implicit_type,special_assumption, pex, pim, p, s, num2str(k));
            f = dir([direc '*.mat']);
            rk = load([direc f(end).name],'X','r');
            [con, ceq] = nlc_imex(rk.X, s, k, pex, pim, p, implicit_type, special_assumption,enforce_positivity);
            ACCURACY = -floor(min(real(-log10(ceq))));
            saveCondition = (ACCURACY <= -14);
            saveCondition = saveCondition && ~any(con > 1e-14);
            if saveCondition; C(i+1,j+1) = rk.r; else C(i+1, j+1) = 0; end
        catch err
            C(i+1, j+1) = 0;
        end
    end
end

fprintf('\n\nPex = %d, Pim = %d, K = %d\n\n', pex, pim,k)
C