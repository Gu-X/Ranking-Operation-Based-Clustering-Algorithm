clear all
clc
close all
load example % load the dataset

a=0.9; % recommended value
M=10;  % recommended value

%% Running Ranking Operation-based Clustering Algorithm
[Centres,IDX]=ROC(data,a,M);
