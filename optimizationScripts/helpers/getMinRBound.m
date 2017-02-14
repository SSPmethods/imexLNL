function r = getMinRBound(pathName, k)


kref = 0; %max(k-1,0);

try
    
    if k ~= kref
        
        kname = ['/K' num2str(kref)];
        pathName = [pathName kname];
        
        files = dir([pathName '/*.mat']);
        
        
        method = files(end).name;
        rk = load([pathName '/' method]);
        
        r = rk.r;
        
    else
        r = 2.05;
    end
catch err
    r = nan;
end

end
