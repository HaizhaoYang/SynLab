% This code tests the robustness of synchrosqueezed curvelet transforms
%
% By Haizhao Yang

close all
clear all

ifclean = 0;
redglobal = 10;
NM = 5;%10;%5
Nx = 512; Ny = 512;
fprintf('sampling rate is %d per unit length in space\n',Nx);
ns = NM*(randn(Nx,Ny)+i*randn(Nx,Ny));

%% Set up parameters for the synchrosqueezed curvelet transform
%set up size of image
%set up size of samples in the space
SPx = 256; SPy = 256;
R_low = 0; R_high = 512/sqrt(2);
%set up the size of samples in [-R_high,R_high]*[-R_high,R_high] in the
%phase domain
NB = 256; %2*NB+1 by 2*NB+1
if NB/2>R_high
    NB = 2*R_high;
end
%set up the smallest size of supports of wave packets
rad = 1;  is_real = 0;  is_unif = 1;  typeNUFFT = 1;  epsl = 1e-2;
N = max(Nx,Ny);

%%set up an image
gdx = [0:Nx-1]/Nx;     gdy = [0:Ny-1]/Ny;    [xo,yo] = ndgrid(gdx,gdy);
fff = zeros(Nx,Ny);
amp = 0.05; F1 = 60; F2 = 60;

xx = xo + amp*sin(2*pi*xo);     yy = yo + amp*sin(2*pi*yo);
f1 = exp(2*pi*i*(+F1*xx +F1*yy));
xx = xo + amp*sin(2*pi*xo);     yy = yo + amp*sin(2*pi*yo);
f2 = exp(2*pi*i*(-F2*xx +F2*yy));

snr = 10*log10((var(f1(:)))/var(ns(:)))

pic = figure;imagesc([0,1],[0,1],real(f1));axis image;colormap (1-gray);xlabel('x_1');ylabel('x_2');axis xy;
pic = figure;imagesc([0,1],[0,1],real(f1+ns));axis image;colormap (1-gray);xlabel('x_1');ylabel('x_2');axis xy;

for band = 0:1
    for tcnt = 1:3
        for scnt = 1:tcnt
            if band
                R_low = 20; R_high = 120;
            else
                R_low = 0; R_high = N/2;
            end
            is_cos = 1;
            t_sc = 1/2 + tcnt/8;
            s_sc = 1/2 + scnt/8;
            t2 = 1;
            
            %% Apply the synchrosqueezed curvelet transform to clean data
            
            if ifclean
                fff = f1;
                
                red = 1;
                [ccc aaa bbb] = gdct_fwd_red2(fff,is_real,[SPx SPy],R_high,R_low,rad,is_cos,t_sc,s_sc,red,t2);
                
                ncl = numel(ccc);
                t1 = size(ccc{1}{1},1);
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
                            temp = ccc{g}{cnt1,cnt2};
                            cc(:,:,cnt_nclp) = temp(:,t2);
                            temp = aaa{g}{cnt1,cnt2};
                            aa(:,:,cnt_nclp) = temp(:,t2);
                            temp = bbb{g}{cnt1,cnt2};
                            bb(:,:,cnt_nclp) = temp(:,t2);
                            cnt_nclp = cnt_nclp + 1;
                        end
                    end
                end
                clear aaa bbb ccc;
                kk1 = zeros(size(cc));
                kk2 = zeros(size(cc));
                gud = find(abs(cc)>epsl);
                kk1(gud) = aa(gud)./cc(gud) / (2*pi*i);
                kk2(gud) = bb(gud)./cc(gud) / (2*pi*i);
                kk1 = (real(kk1));
                kk2 = (real(kk2));
                clear aa bb;
                EXT = 10^10;
                temp = ones(size(cc))*EXT;
                dist = sqrt(kk1.^2+kk2.^2);
                good = find(dist<=R_high & abs(kk1)<N/2 & abs(kk2)<N/2 & dist>=R_low);
                clear dist;
                temp(good) = kk1(good);
                kk1 = temp;
                temp(good) = kk2(good);
                kk2 = temp;
                clear good temp;
                
                WB = N/NB;
                ss_energy = zeros(t1,t2,NB,NB);
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
                            ss_energy(a,b,loc1,loc2) = ss_energy(a,b,loc1,loc2) + abs(tc(gud(g))).^2;
                        end
                    end
                end
                clear kk1 kk2 cc;
                keep = (ss_energy>min(0.1*NM^2,0));
                ss_energy(~(keep==1)) = 0;
                clear keep;
                
                kb1 = sum(ss_energy,4);
                kb1 = permute(kb1,[1 3 2]);
                kb1 = transpose(kb1(:,:,1));
                pic = figure;imagesc([0,1],[-N/2,N/2],kb1);axis square;axis xy;
                xlabel('x1');            ylabel('k1'); colorbar;%colormap (1-gray);
                title('energy distribution of clean data, red=1');
            end
            
            %% Apply the synchrosqueezed curvelet transform to noisy data
            fff = ns;%(f1+ns);
            
            red = 1;
            [ccc aaa bbb] = gdct2_fwd_red2(fff,is_real,[SPx SPy],R_high,R_low,rad,is_cos,t_sc,s_sc,red,t2);
            
            ncl = numel(ccc);
            t1 = size(ccc{1}{1},1);
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
                        temp = ccc{g}{cnt1,cnt2};
                        cc(:,:,cnt_nclp) = temp(:,t2);
                        temp = aaa{g}{cnt1,cnt2};
                        aa(:,:,cnt_nclp) = temp(:,t2);
                        temp = bbb{g}{cnt1,cnt2};
                        bb(:,:,cnt_nclp) = temp(:,t2);
                        cnt_nclp = cnt_nclp + 1;
                    end
                end
            end
            clear aaa bbb ccc;
            kk1 = zeros(size(cc));
            kk2 = zeros(size(cc));
            gud = find(abs(cc)>epsl);
            kk1(gud) = aa(gud)./cc(gud) / (2*pi*i);
            kk2(gud) = bb(gud)./cc(gud) / (2*pi*i);
            kk1 = (real(kk1));
            kk2 = (real(kk2));
            clear aa bb;
            EXT = 10^10;
            temp = ones(size(cc))*EXT;
            dist = sqrt(kk1.^2+kk2.^2);
            good = find(dist<=R_high & abs(kk1)<N/2 & abs(kk2)<N/2 & dist>=R_low);
            clear dist;
            temp(good) = kk1(good);
            kk1 = temp;
            temp(good) = kk2(good);
            kk2 = temp;
            clear good temp;
            
            WB = N/NB;
            ss_energy = zeros(t1,t2,NB,NB);
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
                        ss_energy(a,b,loc1,loc2) = ss_energy(a,b,loc1,loc2) + abs(tc(gud(g))).^2;
                    end
                end
            end
            clear kk1 kk2 cc;
            keep = (ss_energy>min(0.1*NM^2,0));
            ss_energy(~(keep==1)) = 0;
            clear keep;
            
            kb1 = sum(ss_energy,4);
            kb1 = permute(kb1,[1 3 2]);
            kb1 = transpose(kb1(:,:,1));
            %close all;
            pic = figure;imagesc([0,1],[-N/2,N/2],kb1);axis square;axis xy;
            xlabel('x1');            ylabel('k1'); colorbar;
            title('energy distribution of noisy data, red=1');
            
            
            %% Apply the synchrosqueezed curvelet transform to noisy data
            red = redglobal;
            
            [ccc aaa bbb] = gdct2_fwd_red2(fff,is_real,[SPx SPy],R_high,R_low,rad,is_cos,t_sc,s_sc,red,t2);
            
            ncl = numel(ccc);
            t1 = size(ccc{1}{1},1);
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
                        temp = ccc{g}{cnt1,cnt2};
                        cc(:,:,cnt_nclp) = temp(:,t2);
                        temp = aaa{g}{cnt1,cnt2};
                        aa(:,:,cnt_nclp) = temp(:,t2);
                        temp = bbb{g}{cnt1,cnt2};
                        bb(:,:,cnt_nclp) = temp(:,t2);
                        cnt_nclp = cnt_nclp + 1;
                    end
                end
            end
            clear aaa bbb ccc;
            kk1 = zeros(size(cc));
            kk2 = zeros(size(cc));
            gud = find(abs(cc)>epsl);
            kk1(gud) = aa(gud)./cc(gud) / (2*pi*i);
            kk2(gud) = bb(gud)./cc(gud) / (2*pi*i);
            kk1 = (real(kk1));
            kk2 = (real(kk2));
            clear aa bb;
            EXT = 10^10;
            temp = ones(size(cc))*EXT;
            dist = sqrt(kk1.^2+kk2.^2);
            good = find(dist<=R_high & abs(kk1)<N/2 & abs(kk2)<N/2 & dist>=R_low);
            clear dist;
            temp(good) = kk1(good);
            kk1 = temp;
            temp(good) = kk2(good);
            kk2 = temp;
            clear good temp;
            
            WB = N/NB;
            ss_energy = zeros(t1,t2,NB,NB);
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
                        ss_energy(a,b,loc1,loc2) = ss_energy(a,b,loc1,loc2) + abs(tc(gud(g))).^2;
                    end
                end
            end
            clear kk1 kk2 cc;
            keep = (ss_energy>min(0.1*NM^2,0));
            ss_energy(~(keep==1)) = 0;
            clear keep;
            
            kb1 = sum(ss_energy,4);
            kb1 = permute(kb1,[1 3 2]);
            kb1 = transpose(kb1(:,:,1));
            %close all;
            pic = figure;imagesc([0,1],[-N/2,N/2],kb1);axis square;axis xy;
            xlabel('x1');            ylabel('k1'); colorbar;
            head = sprintf('energy distribution of noisy data, red=%d',red);
            title(head);%colormap (1-gray);
            
            %% Apply the synchrosqueezed curvelet transform to noisy data
            
            [ccc aaa bbb] = gdct2_fwd_red2(fff,is_real,[SPx SPy],R_high,R_low,rad,is_cos,t_sc,s_sc,red,t2);
            
            ncl = numel(ccc);
            t1 = size(ccc{1}{1},1);
            nclp = 0;
            for g = 1:ncl
                nclp = nclp + numel(ccc{g});
            end
            t2 = 1;
            cc = zeros(t1,t2,nclp);
            aa = zeros(t1,t2,nclp);
            bb = zeros(t1,t2,nclp);
            cnt_nclp = 1;
            for g=1:ncl
                [szmccc,sznccc] = size(ccc{g});
                for cnt1 = 1:szmccc
                    for cnt2 = 1:sznccc
                        temp = ccc{g}{cnt1,cnt2};
                        cc(:,1,cnt_nclp) = temp(:,t2);
                        temp = aaa{g}{cnt1,cnt2};
                        aa(:,1,cnt_nclp) = temp(:,t2);
                        temp = bbb{g}{cnt1,cnt2};
                        bb(:,1,cnt_nclp) = temp(:,t2);
                        cnt_nclp = cnt_nclp + 1;
                    end
                end
            end
            clear aaa bbb ccc;
            WB = N/NB;
            ss_energy = zeros(t1,t2,NB,NB);
            EXT = 10^10;
            for a=1:t1
                for b=1:t2
                    gud = find(abs(cc(a,b,:))>epsl);
                    temp = abs(cc(a,b,gud));
                    temp = temp(:)';
                    [val ord] = sort(temp,'descend');
                    nb = ceil(length(gud)*0.1);
                    nb = 1;
                    pos = (ord(1:nb));
                    kk1 = real(sum(aa(a,b,gud(pos)))./sum(cc(a,b,gud(pos))) / (2*pi*i));
                    kk2 = real(sum(bb(a,b,gud(pos)))./sum(cc(a,b,gud(pos))) / (2*pi*i));
                    dist = sqrt(kk1.^2+kk2.^2);
                    if ~(dist<=R_high & abs(kk1)<N/2 & abs(kk2)<N/2 & dist>=R_low);
                        kk1 = EXT;
                        kk2 = EXT;
                    end
                    
                    if kk1<EXT
                        tc = cc(a,b,gud(pos));
                        tc = tc(:)';
                        loc1 = floor((kk1+N/2)/WB)+1;
                        loc2 = floor((kk2+N/2)/WB)+1;
                        ss_energy(a,b,loc1,loc2) = ss_energy(a,b,loc1,loc2) + sum(abs(tc).^2);
                    end
                end
            end
            clear aa bb dist good temp kk1 kk2 cc;
            keep = (ss_energy>min(0.1*NM^2,0));
            ss_energy(~(keep==1)) = 0;
            clear keep;
            
            kb1 = sum(ss_energy,4);
            kb1 = permute(kb1,[1 3 2]);
            kb1 = transpose(kb1(:,:,1));
            %close all;
            pic = figure;imagesc([0,1],[-N/2,N/2],kb1);axis square;axis xy;
            xlabel('x1');            ylabel('k1'); colorbar;
            head = sprintf('energy distribution of noisy data, red=%d, largest',red);
            title(head);%colormap (1-gray);
            
            %close all;
        end
    end
end