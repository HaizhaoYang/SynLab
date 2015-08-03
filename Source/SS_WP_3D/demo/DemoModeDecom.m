% This code use 3D synchrosqueezed wave packet transforms to decompose
% wave-like components. It shows how to synchrosqueeze wave packet
% coefficients, how to separate them in the sharpened phase space, and how
% to recover each component.
%
% By Haizhao Yang and Lexing Ying

close all
clear all

if(1)
    %%Set up parameters
    %set up size of image
    Nx = 128*2; Ny = 128*2; Nz = 128*2;
    %set up size of samples in the space
    SPx = 16; SPy = 16; SPz = 16;
    R_low = 0; R_high = 64;
    %set up the size of samples in [-R_high,R_high]*[-R_high,R_high]*[-R_high,R_high] in the
    %phase domain 2*NB(1)+1 by 2*NB(2)+1 by 2*NB(3)+1
    NB = [16 16 16]*2;
    pos = find(NB/2>R_high);
    NB(pos) = 2*R_high;
    %set up the smallest size of supports of wave packets
    rad = 4;
    epsl = 1e-4;
    N = max([Nx,Ny,Nz]);
    
    %% Apply the 3D synchrosqueezed wave packet transforms on uniform data
    if (1)
        fprintf('Test 3D Synchrosqueezed wave packet transforms on uniform data\n');
        %------------------------------------------------------------------
        %set up an image
        %generate the spatial grid
        gdx = [0:Nx-1]/Nx;
        gdy = [0:Ny-1]/Ny;
        gdz = [0:Nz-1]/Nz;
        [xo,yo,zo] = ndgrid(gdx,gdy,gdz);
        fff = zeros(Nx,Ny,Nz);
        if(1)
            amp = 0.05;
            F1 = 20;
            F2 = 20;
            
            xx = xo + amp*sin(2*pi*xo);
            yy = yo + amp*sin(2*pi*yo);
            zz = zo + amp*sin(2*pi*zo);
            f1 = exp(2*pi*i*F1*(xx + yy + zz));
            pic=figure;hold on;
            imagesc([0,1],[0,1],real(f1(:,:,1)));
            colormap gray; title('1st component');xlabel('x_1'); ylabel('x_2');
            colorbar;axis square;axis equal;
            hold off;
            
            xx = xo + amp*sin(2*pi*xo);
            yy = yo + amp*sin(2*pi*yo);
            zz = zo + amp*sin(2*pi*zo);
            f2 = exp(2*pi*i*F2*(-xx + yy +zz));
            pic=figure;hold on;
            imagesc([0,1],[0,1],real(f2(:,:,1)));
            colormap gray; title('2nd component');xlabel('x_1'); ylabel('x_2');
            colorbar;axis square;axis equal;
            hold off;
            
            NM = 0;
            %ADD BOTH REAL AND IMAGINARY NOISE
            ns = NM*(randn(Nx,Ny,Nz)+i*randn(Nx,Ny,Nz));
            fff = (f1+f2+ns);
            pic=figure;hold on;
            imagesc([0,1],[0,1],real(fff(:,:,1)));
            colormap gray;title('Sum of two deformed plane waves');xlabel('x_1'); ylabel('x_2');
            colorbar;axis square;axis equal;
            hold off;
            
            clear xx yy zz ns;
        end
        tic;
        dat = wp3_setup(Nx,Ny,Nz,SPx,SPy,SPz,rad,R_low,R_high);
        is_real = 0;
        %perform 3D synchrosqueeze transform
        %sss is the wave packet coefficients
        sss = wp3_fwd(Nx,Ny,Nz,SPx,SPy,SPz,dat,rad, fff,is_real,R_low,R_high);
        
        %aaa, bbb, and ccc are the wave atom coefficients with the
        %derivatives in b_1, b_2 and b_3
        [aaa,bbb,ccc] = wp3_ext(Nx,Ny,Nz,SPx,SPy,SPz,dat,rad, fff,is_real,R_low,R_high);
        %clear dat;
        ncl = numel(sss);
        ss = zeros(SPx,SPy,SPz,ncl);
        aa = zeros(SPx,SPy,SPz,ncl);
        bb = zeros(SPx,SPy,SPz,ncl);
        cc = zeros(SPx,SPy,SPz,ncl);
        for g=1:ncl
            ss(:,:,:,g) = sss{g};
            aa(:,:,:,g) = aaa{g};
            bb(:,:,:,g) = bbb{g};
            cc(:,:,:,g) = ccc{g};
        end
        clear aaa bbb ccc;
        kk1 = zeros(size(ss));%local wave vector in the x direction
        kk2 = zeros(size(ss));%local wave vector in the y direction
        kk3 = zeros(size(ss));%local wave vector in the z direction
        gud = find(abs(ss)>epsl);
        kk1(gud) = aa(gud)./ss(gud) / (2*pi*i);
        kk2(gud) = bb(gud)./ss(gud) / (2*pi*i);
        kk3(gud) = cc(gud)./ss(gud) / (2*pi*i);
        clear aa bb;
        kk1 = round(real(kk1));
        kk2 = round(real(kk2));
        kk3 = round(real(kk3));
        EXT = 10^10;
        szc = size(ss);
        temp = repmat(EXT,szc);
        dist = sqrt(kk1.^2+kk2.^2+kk3.^2);
        good = find(dist<=R_high & dist>R_low & abs(kk1)<Nx/2 & abs(kk2)<Ny/2 & abs(kk3)<Nz/2);
        temp(good) = kk1(good);
        kk1 = temp;
        temp = repmat(EXT,szc);
        temp(good) = kk2(good);
        kk2 = temp;
        temp = repmat(EXT,szc);
        temp(good) = kk3(good);
        kk3 = temp;
        clear dist;
        
        if(1)
            WB = NB/R_high;
            kb = zeros(2*NB(1)+1,2*NB(1)+1,2*NB(3)+1,SPx,SPy,SPz);
            for a=1:SPx
                for b=1:SPy
                    for c=1:SPz
                        tc = ss(a,b,c,:);
                        tk1 = kk1(a,b,c,:);
                        tk2 = kk2(a,b,c,:);
                        tk3 = kk3(a,b,c,:);
                        tc = tc(:);
                        gud = find(tk1<EXT);
                        tk1 = tk1(gud);
                        tk2 = tk2(gud);
                        tk3 = tk3(gud);
                        for g=1:length(gud)
                            loc1 = round( tk1(g)*WB(1)) + NB(1) + 1;
                            loc2 = round( tk2(g)*WB(2)) + NB(2) + 1;
                            loc3 = round( tk3(g)*WB(3)) + NB(3) + 1;
                            kb(loc1,loc2,loc3,a,b,c) = kb(loc1,loc2,loc3,a,b,c) + abs(tc(gud(g))).^2;
                        end
                    end
                end
            end
            keep = (kb>8*NM^2);
        end
        
        tb = zeros(size(kb));%the ss energy
        tb(keep==1) = kb(keep==1);
        clear keep kb;
        toc
        %------------------------------------------------------------------
        %show the ss energy distribution
        %ss_energy = permute(tb,[4 5 6 1 2 3]);
        figure;imagesc([-R_high,R_high],[-R_high, R_high],squeeze(sum(tb(:,:,:,8,8,8),3)));axis square;colorbar;
        title('ss energy distribution at one sample point');
        temp = abs(fftshift(fftn(fff))/N);
        figure;imagesc([-N/2 N/2],[-N/2 N/2],squeeze(sum(temp,3)));axis square;colorbar;
        title('Fourier spectrum of the image');
        %----------------------------------------
        %separate two modes. In this case, the 1st mode has positive k1 value and the second mode has negative k1 mode
        if(1)
            dd = zeros(size(cc));
            ee = zeros(size(cc));
            for a=1:SPx
                for b=1:SPy
                    for c=1:SPz
                        tk1 = kk1(a,b,c,:);
                        tk2 = kk2(a,b,c,:);
                        tk3 = kk3(a,b,c,:);
                        gud = find(tk1<EXT);
                        tk1 = tk1(gud);
                        tk2 = tk2(gud);
                        tk3 = tk3(gud);
                        for g=1:length(gud)
                            loc1 = round( tk1(g)*WB(1)) + NB(1) + 1;
                            %loc2 = round( tk2(g)*WB(2)) + NB(2) + 1;
                            %loc3 = round( tk3(g)*WB(3)) + NB(3) + 1;
                            if (loc1>NB(1))
                                dd(a,b,c,gud(g)) = ss(a,b,c,gud(g));
                            else
                                ee(a,b,c,gud(g)) = ss(a,b,c,gud(g));
                            end
                        end
                    end
                end
            end
            ddd = cell(size(sss));
            eee = cell(size(sss));
            for g=1:ncl
                ddd{g} = dd(:,:,:,g);
                eee{g} = ee(:,:,:,g);
            end
        end
        
        %----------------------------------------
        %generate the two modes back in image domain and compare with the true value
        if(1)
            g1 = wp3_inv(Nx,Ny,Nz,SPx,SPy,SPz,dat,rad,ddd,is_real,R_low,R_high);
            g2 = wp3_inv(Nx,Ny,Nz,SPx,SPy,SPz,dat,rad,eee,is_real,R_low,R_high);
            e1 = f1-g1;
            e2 = f2-g2;
            figure; imagesc([0,1],[0,1],real(g1(:,:,1))); axis square;colormap gray; colorbar; axis xy; title('recovered 1st component');
            figure; imagesc([0,1],[0,1],real(g2(:,:,1))); axis square;colormap gray; colorbar; axis xy; title('recovered 2nd component');
            fprintf(1, 'ell_2 error of recovered 1st component is %d\n', norm(e1(:))/norm(f1(:)));
            fprintf(1, 'ell_2 error of recovered 2nd component is %d\n', norm(e2(:))/norm(f2(:)));
        end
        
        %----------------------------------------
        %display the (projected) phase space energy distribution after synchrosqueezed transform
        if(1)
            kb1 = tb(:,:,:,:,:,1);
            kb1 = squeeze(sum(kb1,2));
            kb1 = squeeze(sum(kb1,2));
            kb2 = tb(:,:,:,:,:,1);
            kb2 = squeeze(sum(kb2,1));
            kb2 = squeeze(sum(kb2,2));
            kb3 = tb(:,:,:,:,:,1);
            kb3 = squeeze(sum(kb3,1));
            kb3 = squeeze(sum(kb3,1));
            
            g1s = [0:SPx-1]/SPx;
            g2s = [0:SPy-1]/SPy;
            g3s1 = [-R_high:round(1/WB(1)):R_high];
            g3s2 = [-R_high:round(1/WB(2)):R_high];
            g3s3 = [-R_high:round(1/WB(3)):R_high];
            [xx,yy,zz] = meshgrid(g1s,g2s,g3s1);
            figure; %Fixing x3 and stacking k2 k3
            xslice = [0:0.25:1]; yslice = [1]; zslice = [-R_high];
            h = slice(xx,yy,zz,permute(kb1,[2 3 1]),xslice,yslice,zslice);      set(h,'EdgeColor','none');
            xlabel('x1');            ylabel('x2');            zlabel('k1'); colorbar; colormap(1-gray);
            title('energy distribution in phase space, x3 fixed, k2 k3 dimension stacked');
            
            [xx,yy,zz] = meshgrid(g1s,g2s,g3s2);
            figure; %Fixing x3 and stacking k1 k3
            xslice = [1]; yslice = [0:0.25:1];; zslice = [-R_high];
            h = slice(xx,yy,zz,permute(kb2,[2 3 1]),xslice,yslice,zslice);      set(h,'EdgeColor','none');
            xlabel('x1');            ylabel('x2');            zlabel('k2'); colorbar; colormap(1-gray);
            title('energy distribution in phase space, x3 fixed, k1 k3 dimension stacked');
            
            [xx,yy,zz] = meshgrid(g1s,g2s,g3s3);
            figure; %Fixing x3 and stacking k1 k2
            xslice = [1]; yslice = [0:0.25:1]; zslice = [-R_high];
            h = slice(xx,yy,zz,permute(kb3,[2 3 1]),xslice,yslice,zslice);      set(h,'EdgeColor','none');
            xlabel('x1');            ylabel('x2');            zlabel('k3'); colorbar; colormap(1-gray);
            title('energy distribution in phase space, x3 fixed, k1 k2 dimension stacked');
        end
    end
end

