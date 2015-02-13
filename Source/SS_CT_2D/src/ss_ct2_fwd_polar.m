function [kb ccc cc kk1 kk2] = SS_ct_polar(num_direction,fff,SPg,NB,rad,is_real,R_low,R_high,epsl,t_sc, s_sc,is_cos,  wedge_length_coarse)
%Input:
%fff is the image
%SPg(1) is the number of samples in the vertical direction in the image
%SPg(2) is the number of samples in the horizontal direction in the image
%NB(1) is the size of radial grid in [R_low R_high]
%NB(2) is the size of angular grid in [0,pi)
%rad is the smallest size of supports of wave packets
%is_real = 0: fff is complex
%is_real = 1: fff is real
%R_low: lower bound of interested frequency
%R_high: upper bound of interested frequency
%epsl: the threshold for small coefficients of general curvelet transform
%Note: The accuracy of the synchrosqueezed transform is determined by the
%parameter epsl.
%   is_cos      Type of the window function
%                   0: C^infinity window function
%                   1: cosine window function
%               [default set to 0]
%   t_sc    scaling parameter for radius
%               [default set to 1-1/8]
%   s_sc    scaling parameter for angle
%               [default set to 1/2+1/8]
%   wedge_length_coarse
%               length of coarsest wedge 
%               [default set to 5]
%
%Output:
%kb the synchrosqueezed energy distribution in a polar coordinate
%   size(kb) = NB, 
%ccc the general curvelet transform coefficients
%cc a tensor storing the general curvelet transform coefficients
%
%kk1: local wave vector in the first direction
%kk2: local wave vector in the second direction
%
%by Haizhao Yang

if nargin < 12, is_cos = 1; end;
if nargin < 10, t_sc = 1 - 1/8; end;
if nargin < 11, s_sc = 1/2 + 1/8; end;
if nargin < 13, wedge_length_coarse = 8; end;

    %set up size of image
    [Nx Ny] = size(fff);
    N = max(Nx,Ny);
    %set up size of samples in the space
    SPx = SPg(1); SPy = SPg(2);
    
    %ccc is the general curvelet coefficients
    %aaa and bbb are the general curvelet coefficients with the derivatives in b_1 and b_2
        ccc = gdct2_fwd(fff,is_real,1,fff,fff,1,[SPx SPy],R_high,R_low,rad,is_cos, t_sc, s_sc, wedge_length_coarse);
        aaa = gdct2_ext_1(fff,is_real,1,fff,fff,1,[SPx SPy],R_high,R_low,rad,is_cos, t_sc, s_sc, wedge_length_coarse);
        bbb = gdct2_ext_2(fff,is_real,1,fff,fff,1,[SPx SPy],R_high,R_low,rad,is_cos, t_sc, s_sc, wedge_length_coarse);
        
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
        EXT = 10^10;
        szc = size(cc);
        % kk1 local wave vector in the first direction
        % kk2 local wave vector in the second direction
        gud = find(abs(cc)>epsl);
        aa = real( aa(gud)./cc(gud) / (2*pi*i));
        bb = real( bb(gud)./cc(gud) / (2*pi*i));
        dist = sqrt(aa.^2+bb.^2);
        good = find(dist<=R_high & dist>=R_low & abs(aa)<Nx/2 & abs(bb)<Ny/2);
    
        kk1 = repmat(EXT,szc);
        kk1(gud(good)) = aa(good);
        kk2 = kk1;
        kk2(gud(good)) = bb(good);
        
        if(1)
            dr = (R_high-R_low)/(NB(1)-1);
            da = pi/(NB(2)-1)/num_direction;
            %WB = NB/R_high;%2*min(N,2*R_high)/NB;
            kb = zeros(NB(1),NB(2),t1,t2);
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
                        r = sqrt(tk1(g)^2+tk2(g)^2);
                        agl = mod(acos(tk2(g)/r),pi/num_direction);
                        loc1 = round((r-R_low)/dr) + 1;
                        loc2 = round(agl/da) +1;
                        kb(loc1,loc2,a,b) = kb(loc1,loc2,a,b) + abs(tc(gud(g))).^2;
                    end
                end
            end
        end
  

