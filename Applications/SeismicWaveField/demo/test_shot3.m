% This code uses 2D SSCT to extract wave-like components
%
% By Haizhao Yang

close all;
clear all;

if(1)
    N = 512;
    
    if 0 % SS wave packet
        t_sc = 0.7;
        s_sc = 0.7;
        epsl = 0.8*1e-2;
        lemda_k = 17;
        lemda = 3/0.8;
        cont = 5;
        rad = 2;
        thres = 0.001;
    end
    
    if 1 % SS general curvelet
        t_sc = 1-1/8;
        s_sc = 1/2+1/8;
        epsl = 0.5*1e-2;
        lemda_k = 17;%18;
        lemda = 3/0.5;
        cont = 5;
        rad = 2;
        thres = 0.001;
    end
    
    if(1)
        %----------------------------------------
        %setup the image
        fff = zeros(N,N);
        if(1)
            load shot3;
            if(1)
                %tmp = shot3;
                L = size(shot3,2);
                tmp = shot3(1:L,:);
                [szm,szn] = size(tmp);
                
                G1 = linspace(0,1,szm);
                T1 = linspace(0,1,N);
                tmp = csape(G1,permute(tmp,[2 1]));
                tmp = permute(ppval(tmp,T1),[2 1]);
                
                G2 = linspace(0,1,szn);
                T2 = linspace(0,1,N);
                tmp = csape(G2,permute(tmp,[1 2]));
                tmp = permute(ppval(tmp,T2),[1 2]);
                
            end
            K = 80;
            fff = fff - mean(fff(:));
            fff = tmp/(max(abs(tmp(:))));
            
            xx = 0:(1/K):1;
            [yyl yyr] = fdct_wrapping_window(xx);
            xx = [yyl ones(1,N-2*K-2) yyr];
            fff = fff.*(xx'*xx);
            pic = figure; imagesc([0,1],[0,1],real(fff)); colormap gray; title('real data');
            axis square;
            colorbar;
            %saveas(pic,'results/shot3_data.fig');
            size(fff)
            
            NM = 1;
        end
        
        %----------------------------------------
        %perform 2D synchrosqueeze transform
        is_real = 1;
        sz = [N N]/8;
        R_high = N/sqrt(2);
        R_low = 0;
        is_cos = 1;
        is_unif = 1;
        typeNUFFT = 1;
        [xo yo] = ndgrid(0:1/N:(N-1)/N);
        %ccc is the wave atom coefficients
        ccc = gdct2_fwd(fff, is_real,is_unif,xo,yo,typeNUFFT,sz,R_high,R_low,rad,is_cos,t_sc,s_sc);
        %aaa and bbb are the wave atom coefficients with the derivatives in b_1 and b_2
        aaa = gdct2_ext_1(fff, is_real,is_unif,xo,yo,typeNUFFT,sz,R_high,R_low,rad,is_cos,t_sc,s_sc);
        bbb = gdct2_ext_2(fff, is_real,is_unif,xo,yo,typeNUFFT,sz,R_high,R_low,rad,is_cos,t_sc,s_sc);
        
        ncl = numel(ccc);
        [t1 t2] = size(ccc{1}{1});
        NG = t1;
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
        
        EXT = 10^10;
        szc = size(cc);
        gud = find(abs(cc)>epsl);
        aa = real( aa(gud)./cc(gud) / (2*pi*i));
        bb = real( bb(gud)./cc(gud) / (2*pi*i));
        good = find(abs(aa)<=min(R_high,N/2) & abs(bb)<=min(R_high,N/2));
        
        kk1 = repmat(EXT,szc);
        kk1(gud(good)) = aa(good);
        kk2 = kk1;
        kk2(gud(good)) = bb(good);
        
        if (1)
            clear aa aaa bb bbb bad gud;
            nbscales = length(ccc);
            dat = cell(nbscales);
            for cnt = 1:nbscales
                dat{cnt} = cell(size(ccc{cnt}));
            end
            clear ccc;
            NB = t1;
            WB = N/NB;
            kb_cls = cell(t1,t2);
            Num_group = zeros(t1,t2);
            %threshold parameter
            %when thres = 100, almost always 2 groups
            tb = [];
            for a=1:t1
                for b=1:t2
                    M=sparse(NB,NB);
                    tc = cc(a,b,:);
                    tk1 = kk1(a,b,:);
                    tk2 = kk2(a,b,:);
                    tc = tc(:);
                    tk1 = tk1(:);
                    tk2 = tk2(:);
                    for g=1:nclp
                        if tk1(g)~=0 | tk2(g)~=0
                            loc1 = floor((tk1(g)+N/2)/WB)+1;
                            loc2 = floor((tk2(g)+N/2)/WB)+1;
                            if(loc1>=1 && loc1<=NB && loc2>=1 && loc2<=NB)
                                M(loc1,loc2) = M(loc1,loc2) + abs(tc(g)).^2;
                            end
                        end
                    end
                    M(find(M<=thres)) = 0;
                    %clustering at each space point
                    [kb_cls{a,b},cls_ctr,Num_group(a,b)] = angle_cls(1,a,b,M,NB,N,epsl,lemda_k,3);
                    if Num_group(a,b)>0
                        tb = [tb cls_ctr];
                    end
                end
            end
            head = sprintf('Number of clusters at each location: threshold = %f',thres);
            figure;imagesc([0 1],[0 1],Num_group);axis square;colorbar; title(head);
        end
        
        %--------------------------------------------------------------
        %begin clustering
        if (1)
            %cont is a parameter for spectral clustering
            theta = asin(epsl*lemda*N/NB);
            is_periodic = 0;
            [cls_result est_group] = dis_agl_cls(1,tb,cont,theta,t1,NB);
            kb_mode = cell(t1,t2);
            cnt = 0;
            for a = 1:t1
                for b= 1:t2
                    kb_mode{a,b} = sparse(NB,NB);
                    if Num_group(a,b)~=0
                        st = cnt + 1;
                        for cnt_group = 0:Num_group(a,b)-1
                            ctr_pt = cls_result(3:4,st+cnt_group);
                            loc = find(kb_cls{a,b}(3,:)==ctr_pt(1)+NB*(ctr_pt(2)-1));
                            pos1 = kb_cls{a,b}(1,loc);
                            pos2 = kb_cls{a,b}(2,loc);
                            kb_mode{a,b}(pos1+(pos2-1)*NB) = cls_result(5,st+cnt_group);
                        end
                        cnt = cnt + Num_group(a,b);
                    end
                end
            end
            est_group2 = max(max(Num_group));
            pt_per_group = zeros(1,est_group);
            for cnt = 1:est_group
                pt_per_group(cnt) = length(find(cls_result(5,:)==cnt));
            end
            est_group2 = est_group2 + 2;
            large_group = zeros(1,est_group2);
            for cnt = 1:est_group2
                [val large_group(cnt)] = max(pt_per_group);
                pt_per_group(large_group(cnt)) = 0;
            end
            
            mode_ddd = cell(est_group,1);
            mode_g1 = cell(est_group,1);
            for cnt_group = 1:est_group2
                cnt = large_group(cnt_group);
                dd = zeros(size(cc));
                for a = 1:t1
                    for b = 1:t2
                        tc = cc(a,b,:);
                        tk1 = kk1(a,b,:);
                        tk2 = kk2(a,b,:);
                        tk1 = tk1(:);
                        tk2 = tk2(:);
                        for g=1:nclp
                            if tk1(g)~=0 | tk2(g)~= 0
                                loc1 = floor((tk1(g)+N/2)/WB)+1;
                                loc2 = floor((tk2(g)+N/2)/WB)+1;
                                if(loc1>=1 && loc1<=NB && loc2>=1 && loc2<=NB)
                                    if size(kb_mode{a,b},2)>0
                                        if (kb_mode{a,b}(loc1,loc2)==cnt)
                                            dd(a,b,g) = cc(a,b,g);
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
                mode_ddd{cnt_group} = dat;
                cnt_g = 1;
                for g=1:ncl
                    [szmccc,sznccc] = size(dat{g});
                    for cnt1 = 1:szmccc
                        for cnt2 = 1:sznccc
                            mode_ddd{cnt_group}{g}{cnt1,cnt2} = dd(:,:,cnt_g);
                            cnt_g = cnt_g + 1;
                        end
                    end
                end
                mode_g1{cnt_group} = gdct2_inv(mode_ddd{cnt_group},N,N,is_real,R_high,R_low,rad,is_cos,t_sc,s_sc);
                head = sprintf('%d recovered component',cnt_group);
                pic=figure; imagesc([0,1],[0,1],real(mode_g1{cnt_group})); colormap gray; title(head);
                colorbar;
                axis square;
                xlabel('x_1'); ylabel('x_2');
                head = sprintf('results/shot3_recovered_%d.fig',cnt_group);
                %saveas(pic,head);
            end
        end
    end
end



