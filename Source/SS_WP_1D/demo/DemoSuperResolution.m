clear all;
close all;
% This code shows the resolution of the frequency estimation of planewaves,
% i.e., the superresoltuion problem.
%
% By Haizhao Yang

if(1)
    %set up data
    N = 1024;
    x = [0:N-1]/N;
    fff = exp(2*pi*i*x*100.002);
    
    %synchrosqueezed wave packet transform
    eps = 1e-4;
    res = 0.001;
    NG = N/4;
    is_real = 1;
    is_unif = 1;
    typeNUFFT = 1;
    R_high = 150;%N/16;
    R_low = 50;%N/32;
    if is_real
        freq_range = [R_low R_high];
        fff = real(fff);
    else
        freq_range = [-R_high R_high];
    end
    rad = 1.5;
    is_cos = 1;
    t_sc = 1;
    epsl = 1e-4;
    xo = x;
    red = 20;
    
    [T_f coef kk] = ss_wp1_fwd(fff,is_real,is_unif,typeNUFFT,NG,xo,R_high,R_low,rad,is_cos,t_sc,red,epsl,res);
    
    figure;imagesc([0 1],freq_range,real(T_f));colorbar;title('synchrosqueezed energy distribution');axis square;
    axis([0.5,0.50002,99.998,100.008]);
end

if(1)
    %set up data
    N = 1024;
    x = [0:N-1]/N;
    fff = exp(2*pi*i*x*100.00002);
    
    %synchrosqueezed wave packet transform
    eps = 1e-4;
    res = 0.00001;
    NG = N/4;
    is_real = 1;
    is_unif = 1;
    typeNUFFT = 1;
    R_high = 120;%N/16;
    R_low = 80;%N/32;
    if is_real
        freq_range = [R_low R_high];
        fff = real(fff);
    else
        freq_range = [-R_high R_high];
    end
    rad = 1.5;
    is_cos = 1;
    t_sc = 1;
    epsl = 1e-4;
    xo = x;
    red = 20;
    
    [T_f coef kk] = ss_wp1_fwd(fff,is_real,is_unif,typeNUFFT,NG,xo,R_high,R_low,rad,is_cos,t_sc,red,epsl,res);
    
    figure;imagesc([0 1],freq_range,real(T_f));colorbar;title('synchrosqueezed energy distribution');axis square;
    axis([0.5,0.500002,99.99998,100.00008]);
end
