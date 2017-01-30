clear all;
close all;

% This code shows how to tune the parameters of the SST to get better
% estimation of oscillatory instantaneous frequencies. We hae clean and
% noizy examples with good and bad parameters.
%
% By Haizhao Yang

if(1)
    %set up data
    N = 3200;
    x = [0:N-1]/N;
    xx = (150*x*5 + cos(2*pi*x*5)/pi);
    f1 = exp(2*pi*i*xx);
    fff = f1;
    figure;plot(x,real(fff));title('signal');
end

if (1)
    %synchrosqueezed wave packet transform
    eps = 1e-4;
    res = 0.1;
    NG = N/4;
    is_real = 1;
    is_unif = 1;
    typeNUFFT = 1;
    R_high = 800;%N/16;
    R_low = 5*140;%N/32;
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
    red = 10;
    
    [T_f coef kk] = ss_wp1_fwd(fff,is_real,is_unif,typeNUFFT,NG,xo,R_high,R_low,rad,is_cos,t_sc,red,epsl,res);
    
    figure;imagesc([0 1],freq_range,real(T_f));colorbar;title('synchrosqueezed energy distribution');axis square;
end

NM = 1;
ns = NM*randn(1,N);

ff = fff;
fff = fff + ns;

if (1)
    %synchrosqueezed wave packet transform
    eps = 1e-1;
    res = 0.2;
    NG = N/4;
    is_real = 1;
    is_unif = 1;
    typeNUFFT = 1;
    R_high = 800;%N/16;
    R_low = 5*140;%N/32;
    if is_real
        freq_range = [R_low R_high];
        fff = real(fff);
    else
        freq_range = [-R_high R_high];
    end
    rad = 1.5;
    is_cos = 1;
    t_sc = 0.6;
    epsl = 1e-1;
    xo = x;
    red = 20;
    
    [T_f coef kk] = ss_wp1_fwd(fff,is_real,is_unif,typeNUFFT,NG,xo,R_high,R_low,rad,is_cos,t_sc,red,epsl,res);
    
    figure;imagesc([0 1],freq_range,real(T_f));colorbar;title('synchrosqueezed energy distribution');axis square;
end

fff = fff.*exp(-2*pi*i*650*x);

if (1)
    %synchrosqueezed wave packet transform
    eps = 1e-1;
    res = 0.2;
    NG = N/4;
    is_real = 1;
    is_unif = 1;
    typeNUFFT = 1;
    R_high =5*30;%N/16;
    R_low = 5*10;%N/32;
    if is_real
        freq_range = [R_low R_high];
        fff = real(fff);
    else
        freq_range = [-R_high R_high];
    end
    rad = 1.5;
    is_cos = 1;
    t_sc = 0.6;
    epsl = 1e-1;
    xo = x;
    red = 20;
    
    [T_f coef kk] = ss_wp1_fwd(fff,is_real,is_unif,typeNUFFT,NG,xo,R_high,R_low,rad,is_cos,t_sc,red,epsl,res);
    
    figure;imagesc([0 1],freq_range,real(T_f));colorbar;title('synchrosqueezed energy distribution');axis square;
end

fff = ff.*exp(-2*pi*i*650*x)+ns;

if (1)
    %synchrosqueezed wave packet transform
    eps = 1e-1;
    res = 0.2;
    NG = N/4;
    is_real = 1;
    is_unif = 1;
    typeNUFFT = 1;
    R_high =5*30;%N/16;
    R_low = 5*10;%N/32;
    if is_real
        freq_range = [R_low R_high];
        fff = real(fff);
    else
        freq_range = [-R_high R_high];
    end
    rad = 1.5;
    is_cos = 1;
    t_sc = 0.6;
    epsl = 1e-1;
    xo = x;
    red = 20;
    
    [T_f coef kk] = ss_wp1_fwd(fff,is_real,is_unif,typeNUFFT,NG,xo,R_high,R_low,rad,is_cos,t_sc,red,epsl,res);
    
    figure;imagesc([0 1],freq_range,real(T_f));colorbar;title('synchrosqueezed energy distribution');axis square;
end
