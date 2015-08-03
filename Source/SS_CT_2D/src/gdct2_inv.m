function x = gdct2_inv(C, N1, N2, is_real, R_high,R_low, rad, is_cos, t_sc, s_sc, wedge_length_coarse)
% gdct2_inv.m - 2D Inverse Generalized Discrete Curvelet Transform for
% gdct_fwd.m
%
% Inputs
%   C           Cell array of curvelet coefficients.
%               C{j}{l}(k1,k2) is the coefficient at
%                   - scale j: integer, from coarsest to finest scale,
%                   - angle l: integer, starts on the right and
%                   increases counter clockwise,
%                   - position k1,k2: both integers, size doesn't change
%   [N1 N2]     The size of the original image
%
% Optional Inputs
%   is_real     Type of the transform
%                   0: complex-valued curvelets
%                   1: real-valued curvelets
%               [default set to 0]
%   [R_low R_high]         The range of interested spectrum      
%   rad         a parameter to adjust the size of supports of curvelets in
%               the frequency domain, rad <= 2.
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
%               [default set to 8]
%
% Outputs
%   x           N1-by-N2 matrix
%               
%by Haizhao Yang

sz = size(C{1}{1});
XX = zeros(N1,N2);
if nargin < 4, is_real = 0; end;
if nargin < 8, is_cos = 1; end;
if nargin < 9, t_sc = 1 - 1/8; end;
if nargin < 10, s_sc = 1/2 + 1/8; end;
if nargin < 11, wedge_length_coarse = 8; end;
if wedge_length_coarse < 4
    wedge_length_coarse = 4;
    display('wedge_length_coarse is too small. Run with wedge_length_coarse = 4.');
end

% Initialization: data structure
R = sqrt(N1^2+N2^2)/2;
if nargin > 4, R = min(R_high,R); end;
if nargin < 6, R_low = 0; end;
if R_low <= wedge_length_coarse
    R_low = 0;
    display('R_low <= wedge_length_coarse, set R_low = 0.');
end
if nargin < 7, rad = 1; end;
if R<=R_low, error('R_low is too large! R_low at most %f',R); end;
if rad > 2, display('rad too large! rad is set up to be 1'); rad = 1; end;
wedge_length_coarse = wedge_length_coarse*rad;
wedge_length = 2;
wedge_end = [R R+2];
while wedge_end(1) > wedge_length_coarse | length(wedge_end)<4
    Length = rad*wedge_end(1)^t_sc/2;
    wedge_length = [Length wedge_length];
    wedge_end = [(wedge_end(1)-Length) wedge_end];
end
wedge_length = [wedge_end(1) wedge_length];


if mod(N1,2) == 0
pN1 = N1/2;
pos_vec_1 = -pN1:(pN1-1);
else
pN1 = (N1-1)/2;
pos_vec_1 = -pN1:pN1;
end
if mod(N2,2) == 0
pN2 = N2/2;
pos_vec_2 = -pN2:(pN2-1);
else
pN2 = (N2-1)/2;
pos_vec_2 = -pN2:pN2;
end
[Pos_index_1, Pos_index_2] = ndgrid(pos_vec_1,pos_vec_2);

if R_low == 0
nbscales = length(wedge_end);

temp = wedge_length(3:end) +wedge_length(2:end-1);
arc_length = temp.^(s_sc/t_sc);
nbangles = round(pi*wedge_end(3:end)./arc_length/2)*4;
nbangles= [1 1 nbangles];
angles = 2*pi./nbangles;

%NN = ceil(max(wedge_length(end-1)+wedge_length(end-2),sqrt(wedge_end(end-1)^2+wedge_end(end-3)^2-2*wedge_end(end-3)*wedge_end(end-1)*cos(angles(end-1)))));

%Polar coordinate
Pos_radius = sqrt(Pos_index_1.^2+Pos_index_2.^2);
pos = find(Pos_radius>0);
Pos_angle = zeros(N1,N2);
Pos_angle(pos) = acos(Pos_index_1(pos)./Pos_radius(pos));
pos = find(Pos_index_2<0);
Pos_angle(pos) = 2*pi-Pos_angle(pos);
Pos_index_1 = Pos_index_1(:);
Pos_index_2 = Pos_index_2(:);
Pos_radius = Pos_radius(:);
Pos_angle = Pos_angle(:);

level = cell(1,nbscales);%the information of fictitious tiling is nonnecessary
temp = Pos_radius;
for cnt = 1:nbscales
    level{cnt} = struct('polar_r',0,'polar_a',0,'index_1',0,'index_2',0,'X_val',0,'lowpass',0,'hipass',0);
    pos = find(temp<=wedge_end(cnt));
    level{cnt}.polar_r = Pos_radius(pos);
    level{cnt}.polar_a = Pos_angle(pos);
    level{cnt}.index_1 = Pos_index_1(pos);
    level{cnt}.index_2 = Pos_index_2(pos);
    level{cnt}.X_val = zeros(length(pos),1);
    if cnt>1
        if ~is_cos
            level{cnt}.lowpass = fdct_wrapping_window((wedge_end(cnt)-level{cnt}.polar_r)/wedge_length(cnt));
        else
            level{cnt}.lowpass = gdct_cos((wedge_end(cnt)-level{cnt}.polar_r)/wedge_length(cnt));
        end
        level{cnt}.hipass = sqrt(1-(level{cnt}.lowpass).^2);
    end
    temp(pos)= 10^10;
end

tiling = cell(1,nbscales-1);
for cnt = 3:nbscales
    tiling{cnt-1} = cell(2,nbangles(cnt));
    for cnt2 = 1:nbangles(cnt)
        pos = find(level{cnt}.polar_a>=(cnt2-1)*angles(cnt) & level{cnt}.polar_a<cnt2*angles(cnt));
        tiling{cnt-1}{1,cnt2} = struct('polar_a',0,'index_1',0,'index_2',0,'X_val',0,'lowpass',0,'hipass',0);
        tiling{cnt-1}{1,cnt2}.polar_a = level{cnt}.polar_a(pos);
        tiling{cnt-1}{1,cnt2}.index_1 = level{cnt}.index_1(pos);
        tiling{cnt-1}{1,cnt2}.index_2 = level{cnt}.index_2(pos);
        tiling{cnt-1}{1,cnt2}.X_val = zeros(size(pos));
        if ~is_cos
            tiling{cnt-1}{1,cnt2}.lowpass = fdct_wrapping_window((tiling{cnt-1}{1,cnt2}.polar_a-(cnt2-1)*angles(cnt))/angles(cnt));
        else
            tiling{cnt-1}{1,cnt2}.lowpass = gdct_cos((tiling{cnt-1}{1,cnt2}.polar_a-(cnt2-1)*angles(cnt))/angles(cnt));
        end
        tiling{cnt-1}{1,cnt2}.hipass = sqrt(1-(tiling{cnt-1}{1,cnt2}.lowpass).^2);
        tiling{cnt-1}{1,cnt2}.lowpass = (tiling{cnt-1}{1,cnt2}.lowpass).*(level{cnt}.lowpass(pos));
        tiling{cnt-1}{1,cnt2}.hipass = (tiling{cnt-1}{1,cnt2}.hipass).*(level{cnt}.lowpass(pos));
        pos = find(level{cnt-1}.polar_a>=(cnt2-1)*angles(cnt) & level{cnt-1}.polar_a<cnt2*angles(cnt));
        tiling{cnt-1}{2,cnt2} = struct('polar_a',0,'index_1',0,'index_2',0,'X_val',0,'lowpass',0,'hipass',0);
        tiling{cnt-1}{2,cnt2}.polar_a = level{cnt-1}.polar_a(pos);
        tiling{cnt-1}{2,cnt2}.index_1 = level{cnt-1}.index_1(pos);
        tiling{cnt-1}{2,cnt2}.index_2 = level{cnt-1}.index_2(pos);
        tiling{cnt-1}{2,cnt2}.X_val = zeros(size(pos));
        if ~is_cos
            tiling{cnt-1}{2,cnt2}.lowpass = fdct_wrapping_window((tiling{cnt-1}{2,cnt2}.polar_a-(cnt2-1)*angles(cnt))/angles(cnt));
        else
            tiling{cnt-1}{2,cnt2}.lowpass = gdct_cos((tiling{cnt-1}{2,cnt2}.polar_a-(cnt2-1)*angles(cnt))/angles(cnt));
        end
        tiling{cnt-1}{2,cnt2}.hipass = sqrt(1-(tiling{cnt-1}{2,cnt2}.lowpass).^2);
        tiling{cnt-1}{2,cnt2}.lowpass = (tiling{cnt-1}{2,cnt2}.lowpass).*(level{cnt-1}.hipass(pos));
        tiling{cnt-1}{2,cnt2}.hipass = (tiling{cnt-1}{2,cnt2}.hipass).*(level{cnt-1}.hipass(pos));
    end
end
tiling{1} = cell(1,1);
%generate C
%fac = sqrt(sz(1)*sz(2));% sqrt(pi*wedge_end(2)^2);
temp = fftshift(fft2(ifftshift(C{1}{1})))/(sqrt(sz(1)*sz(2)));
x1 = mod(level{1}.index_1,sz(1))+1;
x2 = mod(level{1}.index_2,sz(2))+1;
XX(level{1}.index_1+pN1+1+(level{1}.index_2+pN2+1-1)*N1) = temp(x1+(x2-1)*sz(1));
x1 = mod(level{2}.index_1,sz(1))+1;
x2 = mod(level{2}.index_2,sz(2))+1;
XX(level{2}.index_1+pN1+1+(level{2}.index_2+pN2+1-1)*N1) = temp(x1+(x2-1)*sz(1)).*(level{2}.lowpass);
if 1
    for cnt =3:nbscales
        %fac = sqrt(sz(1)*sz(2));% sqrt(angles(cnt)*(wedge_end(cnt)^2-wedge_end(cnt-2)^2)/4);
        for cnt2 = 1:nbangles(cnt)
            temp = fftshift(fft2(ifftshift(C{cnt-1}{cnt2})))/(sqrt(sz(1)*sz(2)));
            x1 = mod(tiling{cnt-1}{1,cnt2}.index_1+pN1+1,sz(1))+1;
            x2 = mod(tiling{cnt-1}{1,cnt2}.index_2+pN2+1,sz(2))+1;
            tiling{cnt-1}{1,cnt2}.X_val = tiling{cnt-1}{1,cnt2}.X_val + temp(x1+(x2-1)*sz(1)).*(tiling{cnt-1}{1,cnt2}.lowpass);
            x1 = mod(tiling{cnt-1}{2,cnt2}.index_1+pN1+1,sz(1))+1;
            x2 = mod(tiling{cnt-1}{2,cnt2}.index_2+pN2+1,sz(2))+1;
            tiling{cnt-1}{2,cnt2}.X_val = tiling{cnt-1}{2,cnt2}.X_val + temp(x1+(x2-1)*sz(1)).*(tiling{cnt-1}{2,cnt2}.lowpass);
            if cnt2 < nbangles(cnt)
                x1 = mod(tiling{cnt-1}{1,cnt2+1}.index_1+pN1+1,sz(1))+1;
                x2 = mod(tiling{cnt-1}{1,cnt2+1}.index_2+pN2+1,sz(2))+1;
                tiling{cnt-1}{1,cnt2+1}.X_val = tiling{cnt-1}{1,cnt2+1}.X_val+ temp(x1+(x2-1)*sz(1)).*(tiling{cnt-1}{1,cnt2+1}.hipass); 
                x1 = mod(tiling{cnt-1}{2,cnt2+1}.index_1+pN1+1,sz(1))+1;
                x2 = mod(tiling{cnt-1}{2,cnt2+1}.index_2+pN2+1,sz(2))+1;
                tiling{cnt-1}{2,cnt2+1}.X_val = tiling{cnt-1}{2,cnt2+1}.X_val + temp(x1+(x2-1)*sz(1)).*(tiling{cnt-1}{2,cnt2+1}.hipass);
            else
                x1 = mod(tiling{cnt-1}{1,1}.index_1+pN1+1,sz(1))+1;
                x2 = mod(tiling{cnt-1}{1,1}.index_2+pN2+1,sz(2))+1;
                tiling{cnt-1}{1,1}.X_val = tiling{cnt-1}{1,1}.X_val+ temp(x1+(x2-1)*sz(1)).*(tiling{cnt-1}{1,1}.hipass);         
                x1 = mod(tiling{cnt-1}{2,1}.index_1+pN1+1,sz(1))+1;
                x2 = mod(tiling{cnt-1}{2,1}.index_2+pN2+1,sz(2))+1;
                tiling{cnt-1}{2,1}.X_val = tiling{cnt-1}{2,1}.X_val + temp(x1+(x2-1)*sz(1)).*(tiling{cnt-1}{2,1}.hipass);
            end
        end
    end
end
for cnt = 3:nbscales
    for cnt2 = 1:nbangles(cnt)
        pos = find(level{cnt}.polar_a>=(cnt2-1)*angles(cnt) & level{cnt}.polar_a<cnt2*angles(cnt));
        level{cnt}.X_val(pos) = level{cnt}.X_val(pos) + tiling{cnt-1}{1,cnt2}.X_val;
        pos = find(level{cnt-1}.polar_a>=(cnt2-1)*angles(cnt) & level{cnt-1}.polar_a<cnt2*angles(cnt));
        level{cnt-1}.X_val(pos) = level{cnt-1}.X_val(pos) + tiling{cnt-1}{2,cnt2}.X_val;
    end
end
for cnt = 2:nbscales
    XX(level{cnt}.index_1+pN1+1+(level{cnt}.index_2+pN2+1-1)*N1) = XX(level{cnt}.index_1+pN1+1+(level{cnt}.index_2+pN2+1-1)*N1) + level{cnt}.X_val; 
end
%XX = XX(end:-1:1,:);
if is_real
    if mod(N1,2) == 0
        vec1 = 2:N1;
    else
        vec1 = 1:N1;
    end
    if mod(N2,2) == 0
        vec2 = 2:N2;
    else
        vec2 = 1:N2;
    end
    XX(vec1,vec2) = real(XX(vec1,vec2) + rot90(XX(vec1,vec2),2)) + i*imag(XX(vec1,vec2) - rot90(XX(vec1,vec2),2));
    
    if mod(N1,2) == 0
        pos1 = N1/2+1;
    else
        pos1 = (N1+1)/2;
    end
    if mod(N2,2) == 0
        pos2 = N2/2+1;
    else
        pos2 = (N2+1)/2;
    end
    XX(pos1,pos2) = XX(pos1,pos2)/2;
end
x = fftshift(ifft2(ifftshift(XX)))*sqrt(N1*N2);
if is_real, x = real(x); end;

else
    
cut_pos = find(wedge_end<R_low);
cut_num = length(cut_pos);
while cut_num < 2
    Length = rad*wedge_end(1)^t_sc/2;
    wedge_end = [(wedge_end(1)-Length) wedge_end];
    wedge_length = wedge_length(2:end);
    wedge_length = [wedge_end(1) Length wedge_length];
    cut_pos = find(wedge_end<R_low);
    cut_num = length(cut_pos);
end
    
%nbscales = length(wedge_end);
temp = wedge_length(3:end) +wedge_length(2:end-1);
arc_length = temp.^(s_sc/t_sc);
nbangles = round(pi*wedge_end(3:end)./arc_length/2)*4;
nbangles= [1 1 nbangles];

temp = wedge_end(cut_num);
wedge_length = wedge_length(cut_num+2:end);
wedge_length = [R_low-min(2,R_low-temp), min(2,R_low-temp), wedge_end(cut_num+1)-R_low, wedge_length];
wedge_end = wedge_end(cut_num+1:end);
wedge_end = [R_low-min(2,R_low-temp), R_low, wedge_end];

nbangles = nbangles(cut_num+1:end);
nbangles = [1 1 nbangles];
angles = 2*pi./nbangles;

nbscales = length(wedge_end);
temp = cell(1,nbscales-1);
for j = 2:nbscales-1
    temp{j} = C{j-1};
end
C = temp;

%NN = ceil(max(wedge_length(end-1)+wedge_length(end-2),sqrt(wedge_end(end-1)^2+wedge_end(end-3)^2-2*wedge_end(end-3)*wedge_end(end-1)*cos(angles(end-1)))));

%Polar coordinate
Pos_radius = sqrt(Pos_index_1.^2+Pos_index_2.^2);
pos = find(Pos_radius>0);
Pos_angle = zeros(N1,N2);
Pos_angle(pos) = acos(Pos_index_1(pos)./Pos_radius(pos));
pos = find(Pos_index_2<0);
Pos_angle(pos) = 2*pi-Pos_angle(pos);
Pos_index_1 = Pos_index_1(:);
Pos_index_2 = Pos_index_2(:);
Pos_radius = Pos_radius(:);
Pos_angle = Pos_angle(:);

level = cell(1,nbscales);%the information of fictitious tiling is nonnecessary
temp = Pos_radius;
for cnt = 1:nbscales
    level{cnt} = struct('polar_r',0,'polar_a',0,'index_1',0,'index_2',0,'X_val',0,'lowpass',0,'hipass',0);
    pos = find(temp<=wedge_end(cnt));
    level{cnt}.polar_r = Pos_radius(pos);
    level{cnt}.polar_a = Pos_angle(pos);
    level{cnt}.index_1 = Pos_index_1(pos);
    level{cnt}.index_2 = Pos_index_2(pos);
    level{cnt}.X_val = zeros(size(pos));
    if cnt>1
        if ~is_cos
            level{cnt}.lowpass = fdct_wrapping_window((wedge_end(cnt)-level{cnt}.polar_r)/wedge_length(cnt));
        else
            level{cnt}.lowpass = gdct_cos((wedge_end(cnt)-level{cnt}.polar_r)/wedge_length(cnt));
        end
        level{cnt}.hipass = sqrt(1-(level{cnt}.lowpass).^2);
    end
    temp(pos)= 10^10;
end

tiling = cell(1,nbscales-1);
for cnt = 3:nbscales
    tiling{cnt-1} = cell(2,nbangles(cnt));
    for cnt2 = 1:nbangles(cnt)
        pos = find(level{cnt}.polar_a>=(cnt2-1)*angles(cnt) & level{cnt}.polar_a<cnt2*angles(cnt));
        tiling{cnt-1}{1,cnt2} = struct('polar_a',0,'index_1',0,'index_2',0,'X_val',0,'lowpass',0,'hipass',0);
        tiling{cnt-1}{1,cnt2}.polar_a = level{cnt}.polar_a(pos);
        tiling{cnt-1}{1,cnt2}.index_1 = level{cnt}.index_1(pos);
        tiling{cnt-1}{1,cnt2}.index_2 = level{cnt}.index_2(pos);
        tiling{cnt-1}{1,cnt2}.X_val = zeros(size(pos));
        if ~is_cos
            tiling{cnt-1}{1,cnt2}.lowpass = fdct_wrapping_window((tiling{cnt-1}{1,cnt2}.polar_a-(cnt2-1)*angles(cnt))/angles(cnt));
        else
            tiling{cnt-1}{1,cnt2}.lowpass = gdct_cos((tiling{cnt-1}{1,cnt2}.polar_a-(cnt2-1)*angles(cnt))/angles(cnt));
        end
        tiling{cnt-1}{1,cnt2}.hipass = sqrt(1-(tiling{cnt-1}{1,cnt2}.lowpass).^2);
        tiling{cnt-1}{1,cnt2}.lowpass = (tiling{cnt-1}{1,cnt2}.lowpass).*(level{cnt}.lowpass(pos));
        tiling{cnt-1}{1,cnt2}.hipass = (tiling{cnt-1}{1,cnt2}.hipass).*(level{cnt}.lowpass(pos));
        pos = find(level{cnt-1}.polar_a>=(cnt2-1)*angles(cnt) & level{cnt-1}.polar_a<cnt2*angles(cnt));
        tiling{cnt-1}{2,cnt2} = struct('polar_a',0,'index_1',0,'index_2',0,'X_val',0,'lowpass',0,'hipass',0);
        tiling{cnt-1}{2,cnt2}.polar_a = level{cnt-1}.polar_a(pos);
        tiling{cnt-1}{2,cnt2}.index_1 = level{cnt-1}.index_1(pos);
        tiling{cnt-1}{2,cnt2}.index_2 = level{cnt-1}.index_2(pos);
        tiling{cnt-1}{2,cnt2}.X_val = zeros(size(pos));
        if ~is_cos
            tiling{cnt-1}{2,cnt2}.lowpass = fdct_wrapping_window((tiling{cnt-1}{2,cnt2}.polar_a-(cnt2-1)*angles(cnt))/angles(cnt));
        else
            tiling{cnt-1}{2,cnt2}.lowpass = gdct_cos((tiling{cnt-1}{2,cnt2}.polar_a-(cnt2-1)*angles(cnt))/angles(cnt));
        end
        tiling{cnt-1}{2,cnt2}.hipass = sqrt(1-(tiling{cnt-1}{2,cnt2}.lowpass).^2);
        tiling{cnt-1}{2,cnt2}.lowpass = (tiling{cnt-1}{2,cnt2}.lowpass).*(level{cnt-1}.hipass(pos));
        tiling{cnt-1}{2,cnt2}.hipass = (tiling{cnt-1}{2,cnt2}.hipass).*(level{cnt-1}.hipass(pos));
    end
end
tiling{1} = cell(1,1);
if 1
    for cnt =3:nbscales
        %fac = sqrt(sz(1)*sz(2));% sqrt(angles(cnt)*(wedge_end(cnt)^2-wedge_end(cnt-2)^2)/4);
        for cnt2 = 1:nbangles(cnt)
            temp = fftshift(fft2(ifftshift(C{cnt-1}{cnt2})))/(sqrt(sz(1)*sz(2)));
            x1 = mod(tiling{cnt-1}{1,cnt2}.index_1+pN1+1,sz(1))+1;
            x2 = mod(tiling{cnt-1}{1,cnt2}.index_2+pN2+1,sz(2))+1;
            tiling{cnt-1}{1,cnt2}.X_val = tiling{cnt-1}{1,cnt2}.X_val + temp(x1+(x2-1)*sz(1)).*(tiling{cnt-1}{1,cnt2}.lowpass);
            x1 = mod(tiling{cnt-1}{2,cnt2}.index_1+pN1+1,sz(1))+1;
            x2 = mod(tiling{cnt-1}{2,cnt2}.index_2+pN2+1,sz(2))+1;
            tiling{cnt-1}{2,cnt2}.X_val = tiling{cnt-1}{2,cnt2}.X_val + temp(x1+(x2-1)*sz(1)).*(tiling{cnt-1}{2,cnt2}.lowpass);
            if cnt2 < nbangles(cnt)
                x1 = mod(tiling{cnt-1}{1,cnt2+1}.index_1+pN1+1,sz(1))+1;
                x2 = mod(tiling{cnt-1}{1,cnt2+1}.index_2+pN2+1,sz(2))+1;
                tiling{cnt-1}{1,cnt2+1}.X_val = tiling{cnt-1}{1,cnt2+1}.X_val+ temp(x1+(x2-1)*sz(1)).*(tiling{cnt-1}{1,cnt2+1}.hipass); 
                x1 = mod(tiling{cnt-1}{2,cnt2+1}.index_1+pN1+1,sz(1))+1;
                x2 = mod(tiling{cnt-1}{2,cnt2+1}.index_2+pN2+1,sz(2))+1;
                tiling{cnt-1}{2,cnt2+1}.X_val = tiling{cnt-1}{2,cnt2+1}.X_val + temp(x1+(x2-1)*sz(1)).*(tiling{cnt-1}{2,cnt2+1}.hipass);
            else
                x1 = mod(tiling{cnt-1}{1,1}.index_1+pN1+1,sz(1))+1;
                x2 = mod(tiling{cnt-1}{1,1}.index_2+pN2+1,sz(2))+1;
                tiling{cnt-1}{1,1}.X_val = tiling{cnt-1}{1,1}.X_val+ temp(x1+(x2-1)*sz(1)).*(tiling{cnt-1}{1,1}.hipass);         
                x1 = mod(tiling{cnt-1}{2,1}.index_1+pN1+1,sz(1))+1;
                x2 = mod(tiling{cnt-1}{2,1}.index_2+pN2+1,sz(2))+1;
                tiling{cnt-1}{2,1}.X_val = tiling{cnt-1}{2,1}.X_val + temp(x1+(x2-1)*sz(1)).*(tiling{cnt-1}{2,1}.hipass);
            end
        end
    end
end
for cnt = 3:nbscales
    for cnt2 = 1:nbangles(cnt)
        pos = find(level{cnt}.polar_a>=(cnt2-1)*angles(cnt) & level{cnt}.polar_a<cnt2*angles(cnt));
        level{cnt}.X_val(pos) = level{cnt}.X_val(pos) + tiling{cnt-1}{1,cnt2}.X_val;
        pos = find(level{cnt-1}.polar_a>=(cnt2-1)*angles(cnt) & level{cnt-1}.polar_a<cnt2*angles(cnt));
        level{cnt-1}.X_val(pos) = level{cnt-1}.X_val(pos) + tiling{cnt-1}{2,cnt2}.X_val;
    end
end
for cnt = 2:nbscales
    XX(level{cnt}.index_1+pN1+1+(level{cnt}.index_2+pN2+1-1)*N1) = XX(level{cnt}.index_1+pN1+1+(level{cnt}.index_2+pN2+1-1)*N1) + level{cnt}.X_val; 
end
%XX = XX(end:-1:1,:);%!!!!!!!!!!!!!!!!!!!!!!!!!check
if is_real
    if mod(N1,2) == 0
        vec1 = 2:N1;
    else
        vec1 = 1:N1;
    end
    if mod(N2,2) == 0
        vec2 = 2:N2;
    else
        vec2 = 1:N2;
    end
    XX(vec1,vec2) = real(XX(vec1,vec2) + rot90(XX(vec1,vec2),2)) + i*imag(XX(vec1,vec2) - rot90(XX(vec1,vec2),2));
    
    if mod(N1,2) == 0
        pos1 = N1/2+1;
    else
        pos1 = (N1+1)/2;
    end
    if mod(N2,2) == 0
        pos2 = N2/2+1;
    else
        pos2 = (N2+1)/2;
    end
    XX(pos1,pos2) = XX(pos1,pos2)/2;
end
x = fftshift(ifft2(ifftshift(XX)))*sqrt(N1*N2);
if is_real, x = real(x); end;

end
