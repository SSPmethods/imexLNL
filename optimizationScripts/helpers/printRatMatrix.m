function printRatMatrix(A)
%clear all; close all; clc

%mth = 'Method/DIRK/G/Pex2/Pim2/Plin5/S5/K0/method_typeG_r1_1.0000000000000_acc_-16.mat';
%rk = load(mth);

s_r = size(A,1);
s_c = size(A,2);

for i = 1:s_r
    for j = 1:s_c
        
        temp = A(i,j);
        if temp < 1e-15
            temp = 0;
        end
        
        fprintf('%s \t & ',strtrim(rats(temp)));
    end
    fprintf('\n');
end