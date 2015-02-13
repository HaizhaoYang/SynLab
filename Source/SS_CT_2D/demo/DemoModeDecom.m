% this code uses synchrosqueezed curvelet transforms to separate 2D
% wave-like components
%
% By Haizhao Yang

close all
clear all

%% Set up parameters for the synchrosqueezed curvelet transform
%set up size of image
Nx = 512; Ny = 512;
%set up size of samples in the space
SPx = 64; SPy = 64;
R_low = 0; R_high = 512/sqrt(2);
%set up the size of samples in [-R_high,R_high]*[-R_high,R_high] in the
%phase domain
NB = 16; %2*NB+1 by 2*NB+1
if NB/2>R_high
    NB = 2*R_high;
end
%set up the smallest size of supports of wave packets
rad = 1;
N = max(Nx,Ny);
epsl = 1e-2;

%% Apply the synchrosqueezed curvelet transform to data sampled in a uniform grid
if (1)
    fprintf('Test 2D Synchrosqueezed curvelet transforms on uniform data\n');
    %------------------------------------------------------------------
    %set up an image
    %generate the spatial grid
    gdx = [0:Nx-1]/Nx;
    gdy = [0:Ny-1]/Ny;
    [xo,yo] = ndgrid(gdx,gdy);
    fff = zeros(Nx,Ny);
    if(1)
        amp = 0.05;
        F1 = 30;
        F2 = 30;
        %The local wave number is about 85
        
        xx = xo + amp*sin(2*pi*xo);
        yy = yo + amp*sin(2*pi*yo);
        f1 = exp(2*pi*i*(+F1*xx +F1*yy));
        pic=figure;hold on;
        imagesc([0,1],[0,1],real(f1));
        %imagesc([0,0.25],[0,0.25],real(f1(193:256,193:256)));
        colormap gray; title('1st component');xlabel('x_1'); ylabel('x_2');
        colorbar;axis square;axis equal;
        %rectangle('Position',[193/512,193/512,(256-193+1)/512,(256-193+1)/512]);
        hold off;
        
        xx = xo + amp*sin(2*pi*xo);
        yy = yo + amp*sin(2*pi*yo);
        f2 = exp(2*pi*i*(-F2*xx +F2*yy));
        pic=figure;hold on;
        imagesc([0,1],[0,1],real(f2));
        %imagesc([0,0.25],[0,0.25],real(f2(193:256,193:256)));
        colormap gray; title('2nd component');xlabel('x_1'); ylabel('x_2');
        colorbar;axis square;axis equal;
        %rectangle('Position',[193/512,193/512,(256-193+1)/512,(256-193+1)/512]);
        hold off;
        
        NM = 0;
        %ADD BOTH REAL AND IMAGINARY NOISE
        ns = NM*(randn(Nx,Ny)+i*randn(Nx,Ny));
        fff = (f1+f2+ns);
        pic=figure;hold on;
        imagesc([0,1],[0,1],real(fff));
        %imagesc([0,0.25],[0,0.25],real(fff(193:256,193:256)));
        colormap gray;title('Sum of two deformed plane waves');xlabel('x_1'); ylabel('x_2');
        colorbar;axis square;axis equal;
        %rectangle('Position',[193/512,193/512,(256-193+1)/512,(256-193+1)/512]);
        hold off;
        
        is_real = 0;
    end
    
    %----------------------------------------
    %perform 2D synchrosqueeze transform
    is_unif = 1;
    typeNUFFT = 1;
    
    %ccc is the wave atom coefficients
    ccc = gdct2_fwd(fff,is_real,is_unif,xo,yo,typeNUFFT,[SPx SPy],R_high,R_low,rad);
    %aaa and bbb are the wave atom coefficients with the derivatives in b_1 and b_2
    aaa = gdct2_ext_1(fff,is_real,is_unif,xo,yo,typeNUFFT,[SPx SPy],R_high,R_low,rad);
    bbb = gdct2_ext_2(fff,is_real,is_unif,xo,yo,typeNUFFT,[SPx SPy],R_high,R_low,rad);
    
    ncl = numel(ccc);
    [t1 t2] = size(ccc{1}{1});
    nclp = 0;
    for g = 1:ncl
        nclp = nclp + numel(ccc{g});
    end
    cc = zeros(t1,t2,nclp);
    aa = zeros(t1,t2,nclp);
    bb = zeros(t1,t2,nclp);
    cnt_nclp = 1;
    for g=1:ncl
        [szmccc,sznccc] = size(ccc{g});
        for cnt1 = 1:szmccc
            for cnt2 = 1:sznccc
                cc(:,:,cnt_nclp) = ccc{g}{cnt1,cnt2};
                aa(:,:,cnt_nclp) = aaa{g}{cnt1,cnt2};
                bb(:,:,cnt_nclp) = bbb{g}{cnt1,cnt2};
                cnt_nclp = cnt_nclp + 1;
            end
        end
    end
    kk1 = zeros(size(cc));
    kk2 = zeros(size(cc));
    gud = find(abs(cc)>epsl);
    kk1(gud) = aa(gud)./cc(gud) / (2*pi*i);
    kk2(gud) = bb(gud)./cc(gud) / (2*pi*i);
    kk1 = round(real(kk1));
    kk2 = round(real(kk2));
    
    EXT = 10^10;
    temp = ones(size(cc))*EXT;
    dist = sqrt(kk1.^2+kk2.^2);
    good = find(dist<=R_high & abs(kk1)<N/2 & abs(kk2)<N/2 & dist>=R_low);
    temp(good) = kk1(good);
    kk1 = temp;
    temp(good) = kk2(good);
    kk2 = temp;
    
    if(1)
        NB = 64;
        WB = N/NB;
        kb = zeros(t1,t2,NB,NB);
        for a=1:t1
            for b=1:t2
                tc = cc(a,b,:);
                tk1 = kk1(a,b,:);
                tk2 = kk2(a,b,:);
                tc = tc(:);
                gud = find(tk1<EXT);
                tk1 = tk1(gud);
                tk2 = tk2(gud);
                for g=1:length(gud)
                    loc1 = floor((tk1(g)+N/2)/WB)+1;
                    loc2 = floor((tk2(g)+N/2)/WB)+1;
                    kb(a,b,loc1,loc2) = kb(a,b,loc1,loc2) + abs(tc(gud(g))).^2;
                end
            end
        end
        keep = (kb>8*NM^2);
        tb = zeros(size(kb));%the ss energy
        tb(keep==1) = kb(keep==1);
    end
    %%
    %----------------------------------------
    %separate two modes. In this case, the 1st mode has positive k1 value and the second mode has negative k1 mode
    if(1)
        dd = zeros(size(cc));
        ee = zeros(size(cc));
        for a=1:t1
            for b=1:t2
                tk1 = kk1(a,b,:);
                tk2 = kk2(a,b,:);
                gud = find(tk1<EXT);
                tk1 = tk1(gud);
                tk2 = tk2(gud);
                for g=1:length(gud)
                    loc1 = floor((tk1(g)+N/2)/WB)+1;
                    loc2 = floor((tk2(g)+N/2)/WB)+1;
                    if (loc1>NB/2)
                        dd(a,b,gud(g)) = cc(a,b,gud(g));
                    else
                        ee(a,b,gud(g)) = cc(a,b,gud(g));
                    end
                end
            end
        end
        ddd = cell(size(ccc));
        eee = cell(size(ccc));
        cnt_nclp = 1;
        for g=1:ncl
            [szmccc,sznccc] = size(ccc{g});
            for cnt1 = 1:szmccc
                for cnt2 = 1:sznccc
                    ddd{g}{cnt1,cnt2} = dd(:,:,cnt_nclp);
                    eee{g}{cnt1,cnt2} = ee(:,:,cnt_nclp);
                    cnt_nclp = cnt_nclp + 1;
                end
            end
        end
    end
    
    %----------------------------------------
    %generate the two modes back in image domain and compare with the true value
    
    if(1)
        g1 = gdct2_inv(ddd,Nx,Ny,is_real,R_high,R_low, rad);
        g2 = gdct2_inv(eee,Nx,Ny,is_real,R_high,R_low, rad);
        e1 = f1-g1;
        e2 = f2-g2;
        figure; imagesc([0,1],[0,1],real(g1)); colormap gray; colorbar; title('recovered 1st component');
        figure; imagesc([0,1],[0,1],real(g2)); colormap gray; colorbar; title('recovered 2nd component');
        fprintf(1, 'ell_2 error of recovered 1st component is %d\n', norm(e1(:))/norm(f1(:)));
        fprintf(1, 'ell_2 error of recovered 2nd component is %d\n', norm(e2(:))/norm(f2(:)));
    end
    
    %----------------------------------------
    %display the (projected) phase space energy distribution after synchrosqueeze transform
    
    if(1)
        tb = zeros(size(kb));
        tb(keep==1) = kb(keep==1);
        
        kb1 = squeeze(sum(tb,4));
        kb2 = squeeze(sum(tb,3));
        
        g1s = [0:t1-1]/t1;
        g2s = [0:t2-1]/t2;
        g3s = [-N/2:WB:N/2-1];
        [xx,yy,zz] = meshgrid(g1s,g2s,g3s);
        figure;
        xslice = [1]; yslice = [0:0.25:1]; zslice = [-N/2];
        h = slice(xx,yy,zz,permute(kb1,[2 1 3]),xslice,yslice,zslice);      set(h,'EdgeColor','none');
        xlabel('x1');            ylabel('x2');            zlabel('k1'); colorbar; colormap(1-gray);
        title('energy distribution in phase space, k2 dimension stacked');
        
        figure;
        xslice = [0:0.25:1]; yslice = [1]; zslice = [-N/2];
        h = slice(xx,yy,zz,permute(kb2,[2 1 3]),xslice,yslice,zslice);      set(h,'EdgeColor','none');
        xlabel('x1');            ylabel('x2');            zlabel('k2'); colorbar; colormap(1-gray);
        title('energy distribution in phase space, k1 dimension stacked');
    end
end


