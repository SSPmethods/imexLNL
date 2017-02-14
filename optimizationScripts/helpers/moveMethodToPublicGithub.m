% File: moveMethodToPublicGithub.m
% Author: Sidafa Conde
% Email: sconde@umassd.edu
% School: UMass Dartmouth
% Date: 02/02/2017
% Purpose: copy and save optimal methods in the Github folder for the paper, under the SSPMethods organization

clear all; close all; clc

githubPath = '../../imexLNL/';

methodPath = 'Method/DIRK/G';
k = 0;
for pex =2:4
    for pim = pex
        for plin = 2:7
            s = plin + 0;
            kname = num2str(inf);
            
            
            oldName = sprintf('%s/Pex%d/Pim%d/Plin%d/S%d/K0/',...
                methodPath, pex, pim, plin, s)
            
            listOfmethod = dir([oldName '*.mat']);
            
            if ~isempty(listOfmethod)
                bestMethod = sprintf('%s%s',oldName, listOfmethod(end).name);
                
                rk = load(bestMethod);
                
                sigalNaming = sprintf('%ds%dpe%dpi%dplinK%s.mat',...
                    s, pex,pim, plin, kname);
                
                newname = sprintf('%s%s',githubPath, sigalNaming);
                saveGithubMethod(newname, rk.X, rk.r, rk.A, rk.b, rk.At, rk.bt, ...
                rk.s, rk.pex, rk.pim, rk.plin);
            end
        end
    end
end
