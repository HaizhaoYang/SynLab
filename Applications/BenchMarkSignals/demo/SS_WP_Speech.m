clear all;
close all;


%% set up data
load data/speech.txt;
N = length(speech);

x = [0:N-1]/N;
fff = speech(1:N);
fff = fff(:).';
fff = fff/max(fff);

figure;plot(real(fff));

%% synchrosqueezed wave packet transform

% parameters to be tuned
epsl = 1e-6;         % threshold for SST
res =1;             % visualization resolution parameter in frequency
NG = round(N/4);    % number of subsampling points in space
is_real = 1;        % 1: real signals, 0: complex signals
rad = 0.4;          % rad in [0,2] to contral the size of supports of wave packets in the frequency domain
t_sc = 1/2 + 4/8;        % t_sc in [0.5,1], scaling parameter, t_sc = 1 means wavelets, t_sc = 0.5 means wave atoms

% parameters by default
R_high = N/2       % range of interest is [R_low, R_high]
R_low = 0;
is_cos = 1;
is_unif = 1;
typeNUFFT = 1;
red = 32;
sampleRate = 8000;
timeEnd = N/sampleRate;
upBound = round(sampleRate/2);
lowBound = 0;%max(1,round(R_low/res));

% SS wavelet
if 1
    tic;
    [ss_energy coefTensor InsFreq] = ss_wp1_fwd(fff,is_real,is_unif,typeNUFFT,x,NG,R_high,R_low,rad,is_cos,t_sc,red,epsl,res,0);
    toc
    
    %% show results
    % low and high frequency part
    pic = figure;imagesc([0 timeEnd],[0 res*upBound],(ss_energy));
    axis square;xlabel('Time (Second)');ylabel('Frequency (Hz)');axis xy;colormap (1-gray);
    colorbar off;caxis([0 0.5]);
    saveas(pic,'results/Speech_SS_wl.fig');
    set(gca, 'FontSize', 18);
    b=get(gca);
    set(b.XLabel, 'FontSize', 18);set(b.YLabel, 'FontSize', 18);set(b.ZLabel, 'FontSize', 18);set(b.Title, 'FontSize', 18);
    str = 'results/Speech_SS_wl';
    print(gcf, '-depsc2', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    
    % only high frequency part
    [m n] = size(ss_energy);
    pic = figure;imagesc([0 timeEnd],[round(res*upBound/3) res*upBound],(ss_energy(round(m/3):m,:)));
    axis square;xlabel('Time (Second)');ylabel('Frequency (Hz)');axis xy;colormap (1-gray);
    colorbar off;caxis([0 0.1]);
    saveas(pic,'results/Speech_SS_wl_h.fig');
    set(gca, 'FontSize', 18);
    b=get(gca);
    set(b.XLabel, 'FontSize', 18);set(b.YLabel, 'FontSize', 18);set(b.ZLabel, 'FontSize', 18);set(b.Title, 'FontSize', 18);
    str = 'results/Speech_SS_wl_h';
    print(gcf, '-depsc2', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
end

% SS wave packet
if 1
    rad = 0.8;          % rad in [0,2] to contral the size of supports of wave packets in the frequency domain
    t_sc = 1/2 + 2/8;
    tic;
    [ss_energy coefTensor InsFreq] = ss_wp1_fwd(fff,is_real,is_unif,typeNUFFT,x,NG,R_high,R_low,rad,is_cos,t_sc,red,epsl,res,0);
    toc
    
    %% show results
    % low and high frequency part
    pic = figure;imagesc([0 timeEnd],[0 res*upBound],(ss_energy));
    axis square;xlabel('Time (Second)');ylabel('Frequency (Hz)');axis xy;colormap (1-gray);
    colorbar off;caxis([0 0.5]);
    saveas(pic,'results/Speech_SS_wp_2.fig');
    set(gca, 'FontSize', 18);
    b=get(gca);
    set(b.XLabel, 'FontSize', 18);set(b.YLabel, 'FontSize', 18);set(b.ZLabel, 'FontSize', 18);set(b.Title, 'FontSize', 18);
    str = 'results/Speech_SS_wp_2';
    print(gcf, '-depsc2', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    
    % only high frequency part
    [m n] = size(ss_energy);
    pic = figure;imagesc([0 timeEnd],[round(res*upBound/3) res*upBound],(ss_energy(round(m/3):m,:)));
    axis square;xlabel('Time (Second)');ylabel('Frequency (Hz)');axis xy;colormap (1-gray);
    colorbar off;caxis([0 0.1]);
    saveas(pic,'results/Speech_SS_wp_2h.fig');
    set(gca, 'FontSize', 18);
    b=get(gca);
    set(b.XLabel, 'FontSize', 18);set(b.YLabel, 'FontSize', 18);set(b.ZLabel, 'FontSize', 18);set(b.Title, 'FontSize', 18);
    str = 'results/Speech_SS_wp_2h';
    print(gcf, '-depsc2', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
end

% SS wave packet
if 1
    rad = 1.5;          % rad in [0,2] to contral the size of supports of wave packets in the frequency domain
    t_sc = 1/2 + 1/8;
    tic;
    [ss_energy coefTensor InsFreq] = ss_wp1_fwd(fff,is_real,is_unif,typeNUFFT,x,NG,R_high,R_low,rad,is_cos,t_sc,red,epsl,res,0);
    toc
    
    %% show results
    % low and high frequency part
    pic = figure;imagesc([0 timeEnd],[0 res*upBound],(ss_energy));
    axis square;xlabel('Time (Second)');ylabel('Frequency (Hz)');axis xy;colormap (1-gray);
    colorbar off;caxis([0 0.5]);
    saveas(pic,'results/Speech_SS_wp_1.fig');
    set(gca, 'FontSize', 18);
    b=get(gca);
    set(b.XLabel, 'FontSize', 18);set(b.YLabel, 'FontSize', 18);set(b.ZLabel, 'FontSize', 18);set(b.Title, 'FontSize', 18);
    str = 'results/Speech_SS_wp_1';
    print(gcf, '-depsc2', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    
    % only high frequency part
    [m n] = size(ss_energy);
    pic = figure;imagesc([0 timeEnd],[round(res*upBound/3) res*upBound],(ss_energy(round(m/3):m,:)));
    axis square;xlabel('Time (Second)');ylabel('Frequency (Hz)');axis xy;colormap (1-gray);
    colorbar off;caxis([0 0.1]);
    saveas(pic,'results/Speech_SS_wp_1h.fig');
    set(gca, 'FontSize', 18);
    b=get(gca);
    set(b.XLabel, 'FontSize', 18);set(b.YLabel, 'FontSize', 18);set(b.ZLabel, 'FontSize', 18);set(b.Title, 'FontSize', 18);
    str = 'results/Speech_SS_wp_1h';
    print(gcf, '-depsc2', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
end
