clear all;
close all;


%% set up data
load data/seismic_trace.txt;
% data from: http://library.seg.org/doi/abs/10.1190/geo2013-0204.1
% Roberto H. Herrera, Jiajun Han, and Mirko van der Baan (2014). 
% Applications of the synchrosqueezing transform in seismic time-frequency analysis, 79(3), V55-V64.

% Advantages of SS wave packet transform: we can see clear components in
% the high frequency part and the estimates of low frequency compoents are
% still accurate. If SS wavelet is applied, it would mix high frequency
% components together.

N = length(seismic_trace);

x = [0:N-1]/N;
fff = seismic_trace(1:N);
fff = fff(:).';
fff = fff/max(fff);

figure;plot(real(fff));

%% synchrosqueezed wave packet transform

% parameters by default
is_cos = 1;
is_unif = 1;
typeNUFFT = 1;

% parameters to be tuned
epsl = 1e-6;         % threshold for SST
res =0.5;             % visualization resolution parameter in frequency
NG = round(N/2);    % number of subsampling points in space
is_real = 1;        % 1: real signals, 0: complex signals
% rad in [0,2] to contral the size of supports of wave packets in the frequency domain
sampleRate = 500;
timeEnd = N/sampleRate;
upBound = round(sampleRate/2/res);
lowBound = 0;%max(1,round(R_low/res));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters 
t_sc = 1/2;    % t_sc in [0.5,1], scaling parameter, t_sc = 1 means wavelets, t_sc = 0.5 means wave atoms
red = 32;
rad = 1;
R_high = round(N/2*150/250);       % range of interest is [R_low, R_high]
R_low = 0;
[ss_energy coefTensor InsFreq] = ss_wp1_fwd(fff,is_real,is_unif,typeNUFFT,x,NG,R_high,R_low,rad,is_cos,t_sc,red,epsl,res,0);

% show results
pic = figure;imagesc([0 timeEnd],[0 round(res*upBound*150/250)],ss_energy);axis square;xlabel('Time (Second)');ylabel('Frequency (Hz)');axis xy;
title('t=1/2, SS wave atom');colorbar off;caxis([0 0.5]);colormap (1-gray);
saveas(pic,'results/Trace_SS_wp_0.fig');
set(gca, 'FontSize', 18);
b=get(gca);
set(b.XLabel, 'FontSize', 18);set(b.YLabel, 'FontSize', 18);set(b.ZLabel, 'FontSize', 18);set(b.Title, 'FontSize', 18);
str = 'results/Trace_SS_wp_0';
print(gcf, '-depsc2', str);      command = sprintf('epstopdf %s.eps',str);      system(command);

% high frequency part
R_high = round(N/2*100/250);       % range of interest is [R_low, R_high]
R_low = round(N/2*50/250);
[ss_energy coefTensor InsFreq] = ss_wp1_fwd(fff,is_real,is_unif,typeNUFFT,x,NG,R_high,R_low,rad,is_cos,t_sc,red,epsl,res,0);

% show results
pic = figure;imagesc([0 timeEnd],[round(res*upBound*50/250) round(res*upBound*150/250)],ss_energy);axis square;xlabel('Time (Second)');ylabel('Frequency (Hz)');axis xy;
title('t=1/2, SS wave atom');colorbar off;caxis([0 0.5]);colormap (1-gray);
saveas(pic,'results/Trace_SS_wp_0h.fig');
set(gca, 'FontSize', 18);
b=get(gca);
set(b.XLabel, 'FontSize', 18);set(b.YLabel, 'FontSize', 18);set(b.ZLabel, 'FontSize', 18);set(b.Title, 'FontSize', 18);
str = 'results/Trace_SS_wp_0h';
print(gcf, '-depsc2', str);      command = sprintf('epstopdf %s.eps',str);      system(command);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_sc = 1/2 + 1/8;
% parameters
red = 32;
rad = 0.8;
R_high = round(N/2*150/250);       % range of interest is [R_low, R_high]
R_low = 0;
[ss_energy coefTensor InsFreq] = ss_wp1_fwd(fff,is_real,is_unif,typeNUFFT,x,NG,R_high,R_low,rad,is_cos,t_sc,red,epsl,res,0);

% show results
pic = figure;imagesc([0 timeEnd],[0 round(res*upBound*150/250)],ss_energy);axis square;xlabel('Time (Second)');ylabel('Frequency (Hz)');axis xy;
title('t = 1/2+1/8');colorbar off;caxis([0 0.5]);colormap (1-gray);
saveas(pic,'results/Trace_SS_wp_1.fig');
set(gca, 'FontSize', 18);
b=get(gca);
set(b.XLabel, 'FontSize', 18);set(b.YLabel, 'FontSize', 18);set(b.ZLabel, 'FontSize', 18);set(b.Title, 'FontSize', 18);
str = 'results/Trace_SS_wp_1';
print(gcf, '-depsc2', str);      command = sprintf('epstopdf %s.eps',str);      system(command);

% high frequency part
R_high = round(N/2*100/250);       % range of interest is [R_low, R_high]
R_low = round(N/2*50/250);
[ss_energy coefTensor InsFreq] = ss_wp1_fwd(fff,is_real,is_unif,typeNUFFT,x,NG,R_high,R_low,rad,is_cos,t_sc,red,epsl,res,0);

%% show results
pic = figure;imagesc([0 timeEnd],[round(res*upBound*50/250) round(res*upBound*150/250)],ss_energy);axis square;xlabel('Time (Second)');ylabel('Frequency (Hz)');axis xy;
title('t = 1/2+1/8');colorbar off;caxis([0 0.5]);colormap (1-gray);
saveas(pic,'results/Trace_SS_wp_1h.fig');
set(gca, 'FontSize', 18);
b=get(gca);
set(b.XLabel, 'FontSize', 18);set(b.YLabel, 'FontSize', 18);set(b.ZLabel, 'FontSize', 18);set(b.Title, 'FontSize', 18);
str = 'results/Trace_SS_wp_1h';
print(gcf, '-depsc2', str);      command = sprintf('epstopdf %s.eps',str);      system(command);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters 
t_sc = 1/2 + 3/8;        % t_sc in [0.5,1], scaling parameter, t_sc = 1 means wavelets, t_sc = 0.5 means wave atoms
red = 32;
rad = 0.4;
R_high = round(N/2*150/250);       % range of interest is [R_low, R_high]
R_low = 0;
[ss_energy coefTensor InsFreq] = ss_wp1_fwd(fff,is_real,is_unif,typeNUFFT,x,NG,R_high,R_low,rad,is_cos,t_sc,red,epsl,res,0);

% show results
pic = figure;imagesc([0 timeEnd],[0 round(res*upBound*150/250)],ss_energy);axis square;xlabel('Time (Second)');ylabel('Frequency (Hz)');axis xy;
title('t = 1/2+2/8');colorbar off;caxis([0 0.5]);colormap (1-gray);
saveas(pic,'results/Trace_SS_wp_2.fig');
set(gca, 'FontSize', 18);
b=get(gca);
set(b.XLabel, 'FontSize', 18);set(b.YLabel, 'FontSize', 18);set(b.ZLabel, 'FontSize', 18);set(b.Title, 'FontSize', 18);
str = 'results/Trace_SS_wp_2';
print(gcf, '-depsc2', str);      command = sprintf('epstopdf %s.eps',str);      system(command);

% high frequency part
R_high = round(N/2*100/250);       % range of interest is [R_low, R_high]
R_low = round(N/2*50/250);
[ss_energy coefTensor InsFreq] = ss_wp1_fwd(fff,is_real,is_unif,typeNUFFT,x,NG,R_high,R_low,rad,is_cos,t_sc,red,epsl,res,0);

% show results
pic = figure;imagesc([0 timeEnd],[round(res*upBound*50/250) round(res*upBound*150/250)],ss_energy);axis square;xlabel('Time (Second)');ylabel('Frequency (Hz)');axis xy;
title('t = 1/2+2/8');colorbar off;caxis([0 0.5]);colormap (1-gray);
saveas(pic,'results/Trace_SS_wp_2h.fig');
set(gca, 'FontSize', 18);
b=get(gca);
set(b.XLabel, 'FontSize', 18);set(b.YLabel, 'FontSize', 18);set(b.ZLabel, 'FontSize', 18);set(b.Title, 'FontSize', 18);
str = 'results/Trace_SS_wp_2h';
print(gcf, '-depsc2', str);      command = sprintf('epstopdf %s.eps',str);      system(command);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters 
t_sc = 1;        % t_sc in [0.5,1], scaling parameter, t_sc = 1 means wavelets, t_sc = 0.5 means wave atoms
red = 32;
rad = 0.2;
R_high = round(N/2*150/250);       % range of interest is [R_low, R_high]
R_low = 0;
[ss_energy coefTensor InsFreq] = ss_wp1_fwd(fff,is_real,is_unif,typeNUFFT,x,NG,R_high,R_low,rad,is_cos,t_sc,red,epsl,res,0);

% show results
pic = figure;imagesc([0 timeEnd],[0 round(res*upBound*150/250)],ss_energy);axis square;xlabel('Time (Second)');ylabel('Frequency (Hz)');axis xy;
title('t=1, SS wavelet');colorbar off;caxis([0 0.5]);colormap (1-gray);
saveas(pic,'results/Trace_SS_wl.fig');
set(gca, 'FontSize', 18);
b=get(gca);
set(b.XLabel, 'FontSize', 18);set(b.YLabel, 'FontSize', 18);set(b.ZLabel, 'FontSize', 18);set(b.Title, 'FontSize', 18);
str = 'results/Trace_SS_wl';
print(gcf, '-depsc2', str);      command = sprintf('epstopdf %s.eps',str);      system(command);

% high frequency part
R_high = round(N/2*100/250);       % range of interest is [R_low, R_high]
R_low = round(N/2*50/250);
[ss_energy coefTensor InsFreq] = ss_wp1_fwd(fff,is_real,is_unif,typeNUFFT,x,NG,R_high,R_low,rad,is_cos,t_sc,red,epsl,res,0);

% show results
pic = figure;imagesc([0 timeEnd],[round(res*upBound*50/250) round(res*upBound*150/250)],ss_energy);axis square;xlabel('Time (Second)');ylabel('Frequency (Hz)');axis xy;
title('t=1, SS wavelet');colorbar off;caxis([0 0.5]);colormap (1-gray);
saveas(pic,'results/Trace_SS_wl_h.fig');
set(gca, 'FontSize', 18);
b=get(gca);
set(b.XLabel, 'FontSize', 18);set(b.YLabel, 'FontSize', 18);set(b.ZLabel, 'FontSize', 18);set(b.Title, 'FontSize', 18);
str = 'results/Trace_SS_wl_h';
print(gcf, '-depsc2', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
