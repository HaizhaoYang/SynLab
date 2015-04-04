function [T_f coef kk] = ss_wp1_fwd(x,is_real,is_unif,typeNUFFT,xo,NG,R_high,R_low,rad,is_cos,t_sc,red,epsl,h,is_fac,wedge_length_coarse)
% ss_wp_fwd_1D.m - 1D Synchrosqueezed Wave Packet Transform
%
% Inputs
%   x           input signal, a vector of length N
%
% Optional Inputs
%   is_real     Type of the transform
%                   0: complex-valued wave packets
%                   1: real-valued wave packets
%               [default set to 0]
%   is_unif     whether x is sampled on a uniform grid
%               0: No;  1: Yes
%   typeNUFFT   1: NUFFT by Air Force Lab
%               2: USFFT by E. Candes
%               3: NUFFT by L. Greengard and J.-Y. Lee
%               4: Direct non-uniform Fourier Transform
%   xo          non-uniform locations at which x is measured
%   NG          number of subsampled points in time
%   [R_low R_high]         The range of interested spectrum
%   rad         a parameter to adjust the size of supports of the mother wave packet in
%               the frequency domain, rad <= 2.
%   is_cos      Type of the window function
%                   0: C^infinity window function
%                   1: cosine window function
%               [default set to 0]
%   t_sc        scaling parameter for radius
%               [default set to 1-1/4]
%   red         redundancy parameter, red is a positive integer
%               [default set to 1]
%   epsl        threshold for instantaneous frequency estimates
%               [default set to 1e-2]
%   h           frequency band width per pixel in the synchrosqueezed
%               time-frequency representation
%               [default set to 1]
%   is_fac      0: do not increase the magnitude of high frequency wave
%               packet coefficients; 1: increase;
%               [default set to 1, better to visualize high frequency
%               instantaneous frequencies]
%   wedge_length_coarse
%               length of coarsest wedge
%               [default set to 4]
%
% Outputs
%   T_f         1D synchrosqueezed wave packet transform, a matrix with NG columns, each column represent frequency
%               information at a fixed time
%   coef        1D wave packet transform coefficients of x
%   kk          instantaneous frequency estimates from each wave packet
%               coefficient
%
%by Haizhao Yang

N = length(x);
if nargin<2, is_real = 1; end;
if nargin<3, is_unif = 1; end;
if nargin<4, typeNUFFT = 1; end;
if nargin<5, xo = 0:1/N:(N-1)/N; end;
if nargin<6, NG = N; end;
if nargin<7, R_high = round(N/2); end;
if nargin<8, R_low = 0; end;
if nargin<9, rad = 1; end;
if nargin<10, is_cos = 1; end;
if nargin<11, t_sc = 1-1/4; end;
if nargin<12, red = 1; end;
if nargin<13, epsl = 1e-2; end;
if nargin<14, h = 1; end;
if nargin<15, is_fac = 1; end;
if nargin<16, wedge_length_coarse = 4; end;
%if ~is_real, R_low = 0; end;

ccc = wp1_fwd(x,is_real,is_unif,typeNUFFT,xo,NG,R_high,R_low,rad,is_cos, t_sc, red, is_fac, wedge_length_coarse);
aaa = wp1_ext(x,is_real,is_unif,typeNUFFT,xo,NG,R_high,R_low,rad,is_cos, t_sc, red, is_fac, wedge_length_coarse);

if is_real
    fqscale = [R_low R_high];
else
    fqscale = [-R_high R_high];
end
num_grid = ceil((fqscale(2)-fqscale(1))/h);
h = (fqscale(2)-fqscale(1))/(num_grid-1);
%grid = fqscale(1):h:fqscale(2);

EXT = 10^10;
ncl = cell(1,red);
nclp = cell(1,red);
NG = length(ccc{1}{1}{1});
EXT = 10^10;
kk = cell(1,red);
T_f = zeros(num_grid,NG);
coef = cell(1,red);
for cntred = 1:red
    ncl{cntred} = numel(ccc{cntred});
    nclp{cntred} = 0;
    for cnt = 1:ncl{cntred}
        nclp{cntred} = nclp{cntred} + numel(ccc{cntred}{cnt});
    end
    
    if is_real
        cc = zeros(nclp{cntred}/2,NG);
        aa = zeros(nclp{cntred}/2,NG);
        cnt_nclp = 1;
        for cnt = 1:ncl{cntred}
            cc(cnt_nclp,:) = ccc{cntred}{cnt}{1};
            aa(cnt_nclp,:) = aaa{cntred}{cnt}{1};
            cnt_nclp = cnt_nclp + 1;
        end
    else
        cc = zeros(nclp{cntred}-1,NG);
        aa = zeros(nclp{cntred}-1,NG);
        cnt_nclp1 = 1;
        cnt_nclp2 = 1;
        for cnt = 1:ncl{cntred}
            for cnt2 = 1:2
                if cnt2 == 1
                    cc(cnt_nclp1+nclp{cntred}/2-1,:) = ccc{cntred}{cnt}{cnt2};
                    aa(cnt_nclp1+nclp{cntred}/2-1,:) = aaa{cntred}{cnt}{cnt2};
                    cnt_nclp1 = cnt_nclp1 + 1;
                else if cnt ~= 1
                        cc(-cnt_nclp2+nclp{cntred}/2,:) = ccc{cntred}{cnt}{cnt2};
                        aa(-cnt_nclp2+nclp{cntred}/2,:) = aaa{cntred}{cnt}{cnt2};
                        cnt_nclp2 = cnt_nclp2 + 1;
                    end
                end
            end
        end
    end
    
    szc = size(cc);
    gud = find(abs(cc)>epsl);
    aa=real(aa(gud)./(2*pi*i*cc(gud)));
    if is_real
        good = find(aa<=R_high & aa>=R_low);
    else
        good = find(aa<=R_high & aa>=-R_high);
    end
    
    kk{cntred} = repmat(EXT,szc);
    kk{cntred}(gud(good)) = aa(good);
    
    for b = 1:NG
        tk = kk{cntred}(:,b);
        gud = find(tk<EXT);
        tk = tk(gud);
        pos = round((tk-fqscale(1))/h)+1;
        for cnt = 1:length(pos)
            T_f(pos(cnt),b) = T_f(pos(cnt),b) + abs(cc(gud(cnt),b)).^2;
        end
    end
    coef{cntred} = cc;
end

if is_real
    T_f = 4*T_f;
end











