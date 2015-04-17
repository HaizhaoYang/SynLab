% This code generate synthetic canvas image and apply 2D SST to analyze
%
% By Haizhao Yang

clear all;
close all;

% Generate data
N = 1024;
x = 0:1/N:(1-1/N);
[xo yo] = ndgrid(x);
x = xo + 0.01*sin(4*pi*xo) + 0.01*sin(4*pi*yo);
y = yo + 0.01*sin(4*pi*yo) + 0.01*sin(4*pi*xo);
X = real(exp(2*pi*i*70*x)+exp(2*pi*i*70*y));
save 'results/synthetic.mat' X;


% Analyze data
fastAnalysis('synthetic',1,0.8,1e-4,1,[8,5],15,90,8,1,2);


