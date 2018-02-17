clear all;
close all;

% This a demo for the wave packet transform and the synchrosqueezed wave
% packet transform.

if(1)
    %set up data
    N = 1024*2;
    % assume the signal is defined on [0,1]
    x = [0:N-1]/N;
    amp = 0.05;
    F1 = 80;
    F2 = 50;
    % xx is the nonlinear phase function of the first component
    xx = x + amp*sin(2*pi*x);
    f1 = exp(2*pi*i*F1*2*xx);
    % yy is the nonlinear phase function of the second component
    yy = x + amp*cos(2*pi*x);
    f2 = exp(2*pi*i*F2*yy);
    
    % add noise if necessary
    NM = 0;
    ns = NM*randn(1,N);
    % obtain a superposition of two components
    fff = f1 + f2 + ns;
    % visualize test signal
    figure;plot(x,real(fff));title('signal');
end

if (1)
    %synchrosqueezed wave packet transform
    % set up parameters, please see the comments in ss_wp1_fwd for their
    % definitions
    res = 1;
    NG = N/16;
    is_real = 1;
    is_unif = 1;
    typeNUFFT = 1;
    R_high = N/8;
    R_low = 0;
    if is_real
        freq_range = [R_low R_high];
        fff = real(fff);
    else
        freq_range = [-R_high R_high];
    end
    rad = 1.5;
    is_cos = 1;
    t_sc = 1/2+2/8;
    epsl = 1e-2;
    xo = x;
    red = 1; % in this example, we set the redundancy to be 1
    
    % perform SS wave packet transform
    [T_f,coef,kk] = ss_wp1_fwd(fff,is_real,is_unif,typeNUFFT,NG,xo,R_high,R_low,rad,is_cos,t_sc,red,epsl,res);
    % T_f is a matrix storing the time-frequency representation by SS wave
    % packet transform
    % coef is a cell structure storing the wave packet transform
    % coefficients
    % kk is a cell structure storing the instantaneous frequency estimation
    % by the wave packet transform
    
    figure;imagesc([0 1],freq_range,real(T_f));colorbar;title('synchrosqueezed energy distribution');axis square;
    figure;imagesc(real(coef{1}));colorbar;title('wave packet transform');axis square;
end

if (1)
    %synchrosqueezed wave packet transform
    res = 1;
    NG = N/16;
    is_real = 1;
    is_unif = 1;
    typeNUFFT = 1;
    R_high = N/8;
    R_low = 0;
    if is_real
        freq_range = [R_low R_high];
        fff = real(fff);
    else
        freq_range = [-R_high R_high];
    end
    rad = 1.5;
    is_cos = 1;
    t_sc = 1/2+2/8;
    epsl = 1e-2;
    xo = x;
    red = 10; % in this example, we set the redundancy to be 10, please compare the results for different redundancy 
    
    [T_f,coef,kk] = ss_wp1_fwd(fff,is_real,is_unif,typeNUFFT,NG,xo,R_high,R_low,rad,is_cos,t_sc,red,epsl,res);

    figure;imagesc([0 1],freq_range,real(T_f));colorbar;title('synchrosqueezed energy distribution');axis square;
    figure;imagesc(real(coef{1}));colorbar;title('wave packet transform');axis square;
    
    % if you would like to modify the SS representation and do the
    % transform, for example, a simple decomposition into two components
    % since the instantaneous frequencies are well separated
    T = cell(1,2);
    % the first component for the first half of the time domain
    T{1} = T_f;
    T{1}(1:100,:) = 0;
    % the second component for the second half of the time domain
    T{2} = T_f;
    T{2}(101:end,:) = 0;
    [ffr, amplt, num_mode] = ss_wp1_invT(T,coef,kk,is_real,N,R_high,R_low,rad,is_cos,t_sc,res);
    figure;subplot(2,1,1);plot(real(ffr(1,:))); subplot(2,1,2);plot(real(ffr(2,:)));
end


