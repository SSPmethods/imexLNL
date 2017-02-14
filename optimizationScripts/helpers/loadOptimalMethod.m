function rk = loadOptimalMethod(pex, pim, plin, s, k)

implicit_type = 'DIRK';
special_assumption = 'G';

pathName = sprintf('Method/%s/%s/Pex%d/Pim%d/Plin%d/S%d/K%s/',...
    implicit_type,special_assumption, pex, pim, plin, s,num2str(k));
files = dir([pathName '/*.mat']);

rk = [];
if ~isempty(files)
    method = files(end).name;
    rk = load([pathName '/' method]);
end

end
