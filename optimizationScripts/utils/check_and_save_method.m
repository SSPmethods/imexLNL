saveCondition = (ACCURACY <= -14);

saveCondition = saveCondition && ~any(con > 1e-13);

if saveCondition
    
    % check that b = bt
    isEqualWeight = isequalVectorEps(b, bt);
    
    % check that c = ct
    isEqualAbscissa = isequalVectorEps(c, ct);
    
    fprintf('(New) (%s, %s) r1 = %17.15f,  r2 = %17.15f \n\n\n',implicit_type, special_assumption, r, rt);
    pathName = sprintf('Method/%s/%s/Pex%d/Pim%d/Plin%d/S%d',...
        implicit_type,special_assumption, pex, pim, plin, s);
    
    kname = ['/K' num2str(k)];
    pathName = [pathName kname];
    
    if exist(pathName, 'file') ~= 7
        mkdir(pathName);
    end
    
    if ~isinf(k)
        basefile = sprintf('method_type%s_r1_%015.13f_acc_%d',special_assumption, r, ACCURACY);
    else
        basefile = sprintf('method_type%s_r2_%015.13f_acc_%d',special_assumption, rt, ACCURACY);
    end
    filename = sprintf('%s/%s', pathName, basefile);
    list_method = length(dir(sprintf('%s*.mat', filename)));
    
    if (list_method > 0) &&( list_method <= 10)
        filename = sprintf('%s_duplicate_%03d', filename, list_method);
    elseif (list_method == 0)
         filename = sprintf('%s', filename);
    else
        fprintf('%s \n', 'Method already exist. No longer saving duplicates!');
        duplicate_method = 1;
        repeat =10;
    end
    filename = sprintf('%s.mat',filename);
    fprintf('%s \n', filename);
    save(filename,'A','b','c','s','At','bt','ct','r','rt','k','X',...
        'implicit_type', 'special_assumption','pex','pim',...
        'plin','n','isEqualWeight','isEqualAbscissa','v', 'alpha', 'alpha_hat');
end
