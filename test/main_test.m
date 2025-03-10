%%  PreProcessing
clear, close all
clc
%%  Initializing

addpath('../main functions');
addpath('../functions');
warning('off','MATLAB:polyshape:repairedBySimplify')

%%  Main

EKFSLAM('method','pro', 'rmethod', 'A');
UKFSLAM('method','pro', 'rmethod', 'B');