function [k] = which_one_to_load(pathName, K)

try
    if K < 1e-14
        k = 0;
        return
    else
        if abs(K-1) < 1e-14
            km1 = 0;
        else
            km1 = K - 1;
        end
        
        % optimial method at K
        direc_k = dir([pathName 'K' num2str(K) '/*.mat']);
        r_k = load([pathName 'K' num2str(K) '/' direc_k(end).name],'r');
        
        % optimial method at K-0.1
        direc_km1 = dir([pathName 'K' num2str(km1) '/*.mat']);
        r_km1 = load([pathName 'K' num2str(km1) '/' direc_km1(end).name],'r');
        
        if r_km1.r >= r_k.r
            k = km1;
        else
            k = K;
        end
    end
    
catch err
    k = K-1;
end

if k  < 1e-14
    k = 0;
end

end