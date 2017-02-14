function queryMethods

for pex = 2:4
    for plin = pex:6
        pim = pex;
        queryMethodClass(pex, pim, plin);
    end
end

end

function queryMethodClass(pex, pim, plin)

% check if b=bt for pex, pim, s for range of k
global Krange;
Krange = [0 1 5 10 100 1000 inf];
qryPath = 'queryLogs';

qryPath = sprintf('%s/pex%d-pim%d-plin%d.txt',qryPath,pex,pim,plin);

% pex = 2;
% pim = 2;
% plin = 2;
S = 2:10;
B = nan(numel(Krange),numel(S));
C = nan(numel(Krange),numel(S));

for j = 1:numel(S)
    [bCheck, cCheck] = checkKRange(pex, pim, plin, S(j));
   B(:,j) = bCheck'; C(:,j) = cCheck';
end

bInfo = [Krange' B];
cInfo = [Krange' C];

fid = fopen(qryPath,'w');

fprintf(fid, '----- Equal Weights -----\n');
printInfo(fid,bInfo);

fprintf(fid, '\n\n----- Equal Abscissa -----\n');
printInfo(fid,cInfo);

fclose(fid);
end

function printInfo(fid,bInfo)
hdr = '  K/S   2    3    4    5    6    7    8    9    10';
lin = repmat('-',1,length(hdr));
sfmt = [repmat('%4d ',1,10) ' \n'];
fprintf(fid,'%s\n',hdr);
fprintf(fid,'%s\n',lin);
fprintf(fid, sfmt,bInfo(:,:).');
end

function [ weightCheckVector, abscissaCheckVector] = ...
    checkKRange(pex, pim, plin, s)
global Krange;

weightCheckVector = nan(size(Krange));
abscissaCheckVector = nan(size(Krange));

for i = 1:numel(Krange)
    k = Krange(i);
    
    rk = getOptimalMethod(pex, pim, plin,s, k);
    isEqualWeight = nan;
    isEqualAbscissa = nan;
    
    if ~isempty(rk)
        rk = load(rk);
        % check that b = bt
        isEqualWeight = isequalVectorEps(rk.b, rk.bt);
        
        % check that c = ct
        isEqualAbscissa = isequalVectorEps(rk.c, rk.ct);
%     else
%         keyboard
    end
    
    weightCheckVector(i) = isEqualWeight;
    abscissaCheckVector(i) = isEqualAbscissa;
end

end
