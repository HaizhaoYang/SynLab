clear all;
close all;

if(1)
    %set up data
    N = 1024*2;
    x = [0:N-1]/N;
    fff = zeros(1,N);
    amp = 0.05;
    F1 = 80;
    F2 = 50;
    xx = x + amp*sin(2*pi*x);
   % figure;plot(F1*(1+amp*2*pi*cos(2*pi*x)));
    f1 = exp(2*pi*i*F1*2*xx);% + exp(2*pi*i*F1*4*xx);% + exp(2*pi*i*F1*20*xx);
    yy = x + amp*cos(2*pi*x);
   % figure;plot(F2*(1-amp*2*pi*sin(2*pi*x)));
    f2 = exp(2*pi*i*F2*yy);
    
    NM = 0;
    ns = NM*randn(1,N);
    fff = f1 + f2 + ns;
    %fff = real(fff);
    figure;plot(x,real(fff));title('signal');
end

if (1)
    %synchrosqueezed wave packet transform
    eps = 1e-2;
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
    red = 1;
    
    [T_f coef kk] = ss_wp1_fwd(fff,is_real,is_unif,typeNUFFT,NG,xo,R_high,R_low,rad,is_cos,t_sc,red,epsl,res);

    figure;imagesc([0 1],freq_range,real(T_f));colorbar;title('synchrosqueezed energy distribution');axis square;
    figure;imagesc(real(coef{1}));colorbar;title('wave packet transform');axis square;
    %figure;imagesc(kk);colorbar;title('kk');
    %size(T_f)
end

if (1)
    %synchrosqueezed wave packet transform
    eps = 1e-2;
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
    red = 10;
    
    [T_f coef kk] = ss_wp1_fwd(fff,is_real,is_unif,typeNUFFT,NG,xo,R_high,R_low,rad,is_cos,t_sc,red,epsl,res);

    figure;imagesc([0 1],freq_range,real(T_f));colorbar;title('synchrosqueezed energy distribution');axis square;
    figure;imagesc(real(coef{1}));colorbar;title('wave packet transform');axis square;
    %figure;imagesc(kk);colorbar;title('kk');
    %size(T_f)
end


