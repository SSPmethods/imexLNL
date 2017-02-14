%Purpose: first check if there's an optimal method


% the lower order method info
pathNameLowerOrder = sprintf('Method/%s/%s/Pex%d/Pim%d/Plin%d/S%d/K%d/',...
    implicit_type,special_assumption, pex, pim, plin-1, s-1,k);

filesLowerOrder = dir([pathNameLowerOrder '/*.mat']);

if ~isempty(filesLowerOrder)
    method = filesLowerOrder(end).name;
    lowerOrderRK = load([pathNameLowerOrder '/' method]);
    minr_bound = lowerOrderRK.r;
else
    minr_bound = nan;
end

pathName = sprintf('Method/%s/%s/Pex%d/Pim%d/Plin%d/S%d',...
    implicit_type,special_assumption, pex, pim, plin, s);

kname = ['/K' num2str(k)];
pathName = [pathName kname];

files = dir([pathName '/*.mat']);

if ~isempty(files)
    method = files(end).name;
    rk = load([pathName '/' method]);
    
    X = rk.X;
    
    if ~isinf(k)
        r = rk.r;
    else
        r = rk.rt;
    end
    
    
    if isnan(minr_bound)
        minr_bound = r + 1e-3;
    end
    
    % check that b = bt
    isEqualWeight = isequalVectorEps(rk.b,rk.bt);
    
    % check that c = ct
    isEqualAbscissa = isequalVectorEps(rk.c,rk.ct);
    
    if norm(minr_bound - r) <= 1e-5
        minr_bound = r + 0.2;
    end
    
    
    minr = r + minrPusher/mainCounter;

    restart = 1;
    n = length(X);
        
    fprintf(1, 'Loading and trying to optimize old method...\n');
    fprintf(1, '\t%s r = %5.4f ---> minr_bound = %5.4f \n', ...
        repmat('*',1,7),r, minr_bound);
    
    free_optimize = 0;
    
else
    % no method to start from
    % so starting from beginning
    
    restart = 0;
    r = 1e-4;
    free_optimize = 0;
    minr = 0.001;
    
    
    if isnan(minr_bound)
        minr_bound = 0.01;
    end
    
    
    fprintf(1, 'Starting the optimizer fresh...\n');
    fprintf(1, '\t%s r = %5.4f ---> minr_bound = %5.4f \n', ...
        repmat('*',1,7),0, minr_bound);
end
