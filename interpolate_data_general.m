clc; clear; close all;

dataFolder = uigetdir(pwd,'Select CSV folder');
d = 0.0254; % length scale for normalziation if needed

%% Grid settings
gridStruct.xmin = -2*d;
gridStruct.xmax =  10*d;
gridStruct.ymin =  0*d;
gridStruct.ymax =  5*d;
gridStruct.dx   =  0.0005;
gridStruct.dy   =  0.0005;

%% Column definitions
colStruct.x = 7;
colStruct.y = 8;
colStruct.u = 2;
colStruct.v = 3;
%colStruct.alpha = 6;   % remove this line if alpha doesn't exist

%% Options
options.maxFiles = 100000;
options.interpMethod = 'natural';
options.extrapMethod = 'nearest';

results = interpolateFlowFieldFolder( ...
            dataFolder, ...
            '*.csv', ...
            gridStruct, ...
            colStruct, ...
            options);

save('InterpolatedResults_IHC110G050.mat','-struct','results','-v7.3');