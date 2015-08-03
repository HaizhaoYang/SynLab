% This code tests 2D (synchrosqueezed) fwd and inv (general) curvelet transform with arbitrary number of samples
%
% By Haizhao Yang

close all
clear all

%% Set up tests
test = cell(6,1); rslt = cell(6,1);
test{1,1} = struct('num',0,'is_real',0,'is_unif',1,'typeNUFFT',1,'gdcterr',[],'ssterr',[]);
test{2,1} = struct('num',0,'is_real',1,'is_unif',1,'typeNUFFT',1,'gdcterr',[],'ssterr',[]);
test{3,1} = struct('num',10,'is_real',0,'is_unif',1,'typeNUFFT',1,'gdcterr',[],'ssterr',[]);
test{4,1} = struct('num',10,'is_real',1,'is_unif',1,'typeNUFFT',1,'gdcterr',[],'ssterr',[]);
test{5,1} = struct('num',11,'is_real',0,'is_unif',1,'typeNUFFT',1,'gdcterr',[],'ssterr',[]);
test{6,1} = struct('num',11,'is_real',1,'is_unif',1,'typeNUFFT',1,'gdcterr',[],'ssterr',[]);

%% Begin tests
for cnt = 1:6
    para = test{cnt,1};
    num = para.num;
    is_real = para.is_real;
    is_unif = para.is_unif;
    typeNUFFT = para.typeNUFFT;
    
    Nx = 512-num; Ny = 512-num;
    gdx = [0:Nx-1]/Nx;
    gdy = [0:Ny-1]/Ny;
    [xo,yo] = ndgrid(gdx,gdy);
    amp = 0.05;
    F1 = 180;
    F2 = 180;
    
    xx = xo + amp*sin(2*pi*xo);
    yy = yo + amp*sin(2*pi*yo);
    f1u = exp(2*pi*i*(+F1*xx +F1*yy));
    
    xx = xo + amp*sin(2*pi*xo);
    yy = yo + amp*sin(2*pi*yo);
    f2u = exp(2*pi*i*(-F2*xx +F2*yy));
    
    NM = 0;
    ns = NM*(randn(Nx,Ny)+i*randn(Nx,Ny));
    ffu = (f1u+f2u+ns);
    
    if ~is_unif
        %Add perturbation
        xo = xo + randn(Nx,Ny)/Ny/16;
        yo = yo + randn(Nx,Ny)/Ny/16;
        
        xx = xo + amp*sin(2*pi*xo);
        yy = yo + amp*sin(2*pi*yo);
        f1 = exp(2*pi*i*(+F1*xx +F1*yy));
        
        xx = xo + amp*sin(2*pi*xo);
        yy = yo + amp*sin(2*pi*yo);
        f2 = exp(2*pi*i*(-F2*xx +F2*yy));
        
        ns = NM*(randn(Nx,Ny)+i*randn(Nx,Ny));
        fff = (f1+f2+ns);
    else
        fff = ffu;
    end
    
    if is_real
        ffi = real(fff);
        ffc = real(ffu);
    else
        ffi = fff;
        ffc = ffu;
    end
    xi = xo;
    yi = yo;
    SPx = 64; SPy = 64;
    R_low = 20; R_high = sqrt(Nx^2+Ny^2)/2;
    NB = 16;
    if NB/2>R_high
        NB = 2*R_high;
    end
    N = max(Nx,Ny);
    epsl = 1e-3;
    rad = 1;
    red = 4;
    t_sc = 0.8;
    s_sc = 0.6;
    coefCell = gdct2_fwd_red(ffi,is_real,[SPx SPy],R_high,R_low,rad,1,0.8,0.6,red);
    norm1 = norm(ffi,2);
    norm2 = 0;
    norm3 = sqrt(sum(sum(abs(ffi).^2)));
    L1 = length(coefCell);
    for cnt1 = 1:L1
        L2 = length(coefCell{cnt1});
        for cnt2 = 1:L2
            norm2 = norm2 + sum(sum(abs(coefCell{cnt1}{cnt2}).^2));
        end
    end
    norm2 = sqrt(norm2);
    [norm1 norm2 norm3]
    ffr = gdct2_inv_red(coefCell,Nx,Ny,is_real,R_high,R_low,rad,1,0.8,0.6,red);
    para.gdcterr = [max(abs(ffc(:)-ffr(:))) norm(ffc-ffr)/norm(ffc)];
    figure;subplot(1,2,1);imagesc(real(ffc-ffr));axis image;colorbar;colormap (1-gray); title('Difference');
    subplot(1,2,2);imagesc(imag(ffc-ffr));axis image;colorbar;colormap (1-gray);
    [ss_energy coefCell coefTensor LocWavVecx LocWavVecy] = ss_ct2_fwd(ffi,[SPx SPy],R_low,R_high,NB,...
        rad,is_real,is_unif,typeNUFFT,xi,yi,epsl,t_sc,s_sc,red);
    ffr = ss_ct2_inv(ss_energy,coefCell,coefTensor,LocWavVecx,LocWavVecy,size(ffc),R_low,R_high,...
        NB,rad,is_real,t_sc,s_sc,red);
    para.ssterr = [max(abs(ffc(:)-ffr(:))) norm(ffc-ffr)/norm(ffc)];
    figure;subplot(1,2,1);imagesc(real(ffc-ffr));axis image;colorbar;colormap (1-gray); title('Difference');
    subplot(1,2,2);imagesc(imag(ffc-ffr));axis image;colorbar;colormap (1-gray);
    rslt{cnt,1} = para;
end

for cnt = 1:6
    rslt{cnt,1}
end




