function img = ss_ct2_inv(ss_energy,coefCell,coefTensor,LocWavVecx,LocWavVecy,imgSize,R_low,R_high,NB,rad,is_real,t_sc,s_sc,red,is_cos,wedge_length_coarse)
% ss_ct2_inv.m - 2D Inverse Synchrosqueezed Curvelet Transform
%
% Inputs:
%   ss_energy        synchrosqueezed energy distribution
%   coefCell         a cell structure storing general curvelet transform coefficients
%   coefTensor       a tensor storing the general curvelet transform coefficients
%   LocWavVecx       local wave vector in the first direction
%   LocWavVecy       local wave vector in the second direction
%   imgSize          size of the original image
%
% Optional inputs:
%   R_low            lower bound of interested frequency, not necessary in this code
%   R_high           upper bound of interested frequency
%   NB               (2*NB+1) by (2*NB+1) is the size of samples in the phase in [-R_high,R_high]*[-R_high,R_high]
%   rad              a parameter to adjust the size of supports of curvelets in
%                    the frequency domain, rad <= 2.
%   is_real          0: img is complex
%                    1: img is real
%   t_sc             scaling parameter for radius
%                    [default set to 1-1/8]
%   s_sc             scaling parameter for angle
%                    [default set to 1/2+1/8]
%   red         a parameter for redundancy
%               [ default set to 1]
%   is_cos           Type of the window function
%                    0: C^infinity window function
%                    1: cosine window function
%                    [default set to 0]
%   wedge_length_coarse
%                    length of coarsest wedge 
%                    [default set to 5]
%
% Output:
%   img              given image of size imgSize(1) by imgSize(2)
%
%by Haizhao Yang

% TODO: insert non-uniform transform

if nargin < 6, error('Not enough inputs!'); end;
Nx = imgSize(1);
Ny = imgSize(2);
N = max(Nx,Ny);
%if nargin < 7, R_low = 0; end;
if nargin < 8, R_high = ceil(N/2); end;
if nargin < 9, NB = ceil(R_high/8); end;
if nargin < 10, rad = 1; end;
if nargin < 11, is_real = 0; end;
if nargin < 12, t_sc = 1 - 1/8; end;
if nargin < 13, s_sc = 1/2 + 1/8; end;
if nargin < 14, red = 1; end;
if nargin < 15, is_cos = 1; end;
if nargin < 16, wedge_length_coarse = 8; end;

dd = zeros(size(coefTensor));
[SPx SPy] = size(coefCell{1}{1});
EXT = 10^10;
ncl = numel(coefCell);
WB = NB/R_high;
for a=1:SPx 
    for b=1:SPy
        tc = coefTensor(a,b,:);
        tk1 = LocWavVecx(a,b,:);
        tk2 = LocWavVecy(a,b,:); 
        tc = tc(:);
        gud = find(tk1<EXT);
        tk1 = tk1(gud);
        tk2 = tk2(gud);
        for g=1:length(gud)  
            loc1 = round(tk1(g)*WB)+NB+1;
            loc2 = round(tk2(g)*WB)+NB+1;
            if( ss_energy(loc1,loc2,a,b)>0)
                 dd(a,b,gud(g)) = coefTensor(a,b,gud(g));
            end
        end
    end
                    
end
ddd = cell(size(coefCell));
cnt_nclp = 1;
for g=1:ncl
    [szmccc,sznccc] = size(coefCell{g});
    for cnt1 = 1:szmccc
       for cnt2 = 1:sznccc
           ddd{g}{cnt1,cnt2} = dd(:,:,cnt_nclp);
           cnt_nclp = cnt_nclp + 1;
       end
    end    
end

R_low = 0;
img = gdct2_inv_red(ddd,Nx,Ny,is_real,R_high*sqrt(2),R_low,rad,is_cos,t_sc,s_sc,red,wedge_length_coarse);