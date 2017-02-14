function saveGithubMethod(name, X, r, A, b, At, bt, s, pex, pim, plin)
% File: saveGithubMethod.m
% Author: Sidafa Conde
% Email: sconde@umassd.edu
% School: UMass Dartmouth
% Date: 02/02/2017
% Purpose: save the method variables to file = name

save(name, 'X', 'r', 'A', 'b', 'At', 'bt', 's', 'pex', 'pim', 'plin');
end
