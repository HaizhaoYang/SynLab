function [mode, amplt, num_mode] = ss_wp_inv_1D(cluster, coef, kk, is_real, N, R_high, R_low, rad, is_cos, t_sc, h, is_fac, wedge_length_coarse)
% ss_wp_inv_1D.m - 1D Inverse Synchrosqueezed Wave Packet Transform
%
% Inputs
%   cluster     a cell structure, cluster{i} contains the positions of
%               the nonzero synchrosqueezed energy of the ith mode
%   coef        1D wave packet transform coefficients of the original signal
%   kk          instantaneous frequency estimates from each wave packet
%               coefficient
%
% Optional Inputs
%   is_real     Type of the transform
%                   0: complex-valued wave packets
%                   1: real-valued wave packets
%               [default set to 0]
%   N           the length of the original signal
%   [R_low R_high]         The range of interested spectrum
%   rad         a parameter to adjust the size of supports of the mother wave packet in
%               the frequency domain, rad <= 2.
%   is_cos      Type of the window function
%                   0: C^infinity window function
%                   1: cosine window function
%               [default set to 0]
%   t_sc        scaling parameter for radius
%               [default set to 1-1/4]
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
%   mode        output signals, a matrix of size #signals by N
%   amplt       amplitude estimates of output signals, a matrix of size #signals by N
%   num_mode    estimated number of modes in cluster
%
%by Haizhao Yang

if nargin<4, is_real = 1; end;
if nargin<5, N = size(coef{1},2); end;
if nargin<6, R_high = round(N/2); end;
if nargin<7, R_low = 0; end;
if nargin<8, rad = 1; end;
if nargin<9, is_cos = 1; end;
if nargin<10, t_sc = 1-1/4; end;
if nargin<11, h = 1; end;
if nargin<12, is_fac = 1; end;
if nargin<13, wedge_length_coarse = 4; end;
if ~is_real, R_low = 0; end;

if is_real
    fqscale = [R_low R_high];
else
    fqscale = [-R_high R_high];
end
num_grid = ceil((fqscale(2)-fqscale(1))/h);
h = (fqscale(2)-fqscale(1))/(num_grid-1);
%grid = fqscale(1):h:fqscale(2);

red = length(coef);
ncl = cell(1,red);
nclp = cell(1,red);
for cntred = 1:red
    if is_real
        [sk NG] = size(coef{cntred});
        nclp{cntred} = sk*2;
    else
        [sk NG] = size(coef{cntred});
        nclp{cntred} = sk + 1;
    end
    ncl{cntred} = nclp{cntred}/2;
end

EXT = 10^10;
num_mode = length(cluster);
mode = zeros(num_mode,N);
amplt = zeros(num_mode,N);
for cnt_mode = 1:num_mode
    ccc = cell(1,red);
    for cntred = 1:red
        cc = zeros(sk,NG);
        for b = 1:NG
            tk = kk{cntred}(:,b);
            gud = find(tk<EXT);
            tk = tk(gud);
            pos = round((tk-fqscale(1))/h)+1;
            for cnt_pos = 1:length(pos)
                if find(cluster{cnt_mode}==pos(cnt_pos)+(b-1)*NG)
                    cc(gud(cnt_pos),b) = coef{cntred}(gud(cnt_pos),b);
                end
            end
        end
        ccc{cntred} = cell(1,ncl{cntred});
        if is_real
            for cnt = 1:ncl{cntred}
                ccc{cntred}{cnt}{1} = cc(cnt,:);
                ccc{cntred}{cnt}{2} = zeros(1,NG);
            end
        else
            ccc{cntred}{1}{2} = zeros(1,NG);
            cnt_nclp1 = 1;
            cnt_nclp2 = 1;
            for cnt = 1:ncl{cntred}
                for cnt2 = 1:2
                    if cnt2 == 1
                        ccc{cntred}{cnt}{cnt2} = cc(cnt_nclp1+nclp{cntred}/2-1,:);
                        cnt_nclp1 = cnt_nclp1 + 1;
                    else if cnt ~= 1
                            ccc{cntred}{cnt}{cnt2} = cc(-cnt_nclp2+nclp{cntred}/2,:);
                            cnt_nclp2 = cnt_nclp2 + 1;
                        end
                    end
                end
            end
        end
    end
    [mode(cnt_mode,:) amplt(cnt_mode,:)] = wp1_inv(ccc, N, is_real, R_high,R_low, rad, is_cos, t_sc, is_fac, wedge_length_coarse);
end
