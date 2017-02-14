clear all; close all; clc

addpath('order_cond/');
addpath('utils/');

fid = fopen('allMethodChecks.md','w');

enforce_positivity = 1;

%get the list of markdown files

md_files = dir('DIRK*.md');
%md_files = dir('DIRK-G-pex2-pim2-plin2-s6.md');

for m = 1:length(md_files)
    
    mdFile = md_files(m).name;
    
    % read the file
    filetext = fileread(mdFile);
    
    % search for the line of code that includes 'mat'
    % each line is separated by a newline ('\n')
    expr = '[^\n]*mat[^\n]*';
    fileread_info = regexp(filetext,expr,'match');
    
%     if ~length(fileread_info)
%         delete(mdFile);
%     end
    
    for i = 1:length(fileread_info)
        
        try
            md_mtd = fileread_info(i);
            
            LogicalStr = {'false', 'true'};
            
            method_filename_cat = md_mtd{1};
            
            load(method_filename_cat)
            
            [A, b, c, At, bt, ct, r, rt] = unpack_imex(X, s, k, implicit_type, special_assumption);
            [v, alpha, alpha_hat]  = Butcher2ShuOsher(A, At, b', bt', r, k);
            [con, ceq] = nlc_imex(X, s, k,pex, pim, plin, implicit_type, special_assumption, enforce_positivity);
            ACCURACY = -floor(min(real(-log10(ceq))));
            saveCondition = (ACCURACY <= -14);
            
            saveCondition = saveCondition && ~any(con > 1e-13);
            
            
            fprintf(fid, '%s \t Accurary: %d \t SSP: %s\n',method_filename_cat, ACCURACY,LogicalStr{(~any(con > 1e-13)+1)});
        catch err
            continue
        end
    end
    
end

fclose(fid);