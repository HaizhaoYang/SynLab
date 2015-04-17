function [x amplt] = wp1_inv(C, N, is_real, R_high,R_low, rad, is_cos, t_sc, is_fac, wedge_length_coarse)
% wp1_inv.m - 1D Inverse Wave Packet Transform
%
% Inputs
%   C           Coefficients of 1D wave packet transform
%   N          The size of the original signal
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
%               [default set to 1-1/4]
%   is_fac      0: do not increase the magnitude of high frequency wave
%               packet coefficients; 1: increase;
%               [default set to 1, better to visualize high frequency
%               instantaneous frequencies]
%   wedge_length_coarse
%               length of coarsest wedge
%               [default set to 4]
%
% Outputs
%   x           recovered signal
%   amplt       recovered amplitude, if x is a single wave-like component
%
%by Haizhao Yang

NG = length(C{1}{1}{1});
XX = zeros(1,N);
if nargin < 3, is_real = 0; end;
if nargin < 7, is_cos = 1; end;
if nargin < 8, t_sc = 1 - 1/4; end;
if nargin < 9, is_fac = 1; end;
if nargin < 10, wedge_length_coarse = 4; end;
if wedge_length_coarse < 4
    wedge_length_coarse = 4;
    display('wedge_length_coarse is too small. Run with wedge_length_coarse = 4.');
end

% Initialization: data structure
R = (N-mod(N,2))/2;
if nargin > 3, R = min(R_high,R); end;
if nargin < 5, R_low = 0; end;
if ~is_real, R_low = 0; end;
if R_low <= wedge_length_coarse
    R_low = 0;
    display('R_low <= wedge_length_coarse, set R_low = 0.');
end
if nargin < 6, rad = 1; end;
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


if mod(N,2) == 0
    pN = N/2;
    pos_vec = -pN:(pN-1);
else
    pN = (N-1)/2;
    pos_vec = -pN:pN;
end
Pos_index = pos_vec;

if R_low == 0
    nbscales = length(wedge_end);
    
    %temp = wedge_length(3:end) +wedge_length(2:end-1);
    
    %Polar coordinate
    Pos_radius = abs(Pos_index);
    
    level = cell(1,nbscales);%the information of fictitious tiling is nonnecessary
    %temp = Pos_radius;
    for cnt = 1:nbscales
        level{cnt} = struct('polar_r',0,'index',0,'X_val',0,'lowpass',0,'hipass',0);
        if cnt == 1
            temp = floor(wedge_end(cnt))+1;
            %st = (-temp+2+(N-mod(N,2))/2);
            %ed = (temp-1+(N+2-mod(N,2))/2);
            st = (-temp+1+(N-mod(N,2))/2);
            ed = (temp-2+(N+2-mod(N,2))/2);
            pos = st:ed;%find(temp<=wedge_end(cnt));
        else %if cnt<nbscales-1
            stOld = st;
            edOld = ed;
            temp = floor(wedge_end(cnt))+1;
            %st = max((-temp+2+(N-mod(N,2))/2),1);
            %ed = min((temp-1+(N+2-mod(N,2))/2),N);
            st = max((-temp+1+(N-mod(N,2))/2),1);
            ed = min((temp-2+(N+2-mod(N,2))/2),N);
            pos = [st:(stOld-1) (edOld+1):ed];%find(temp<=wedge_end(cnt));
        end
        level{cnt}.polar_r = Pos_radius(pos);
        level{cnt}.index = Pos_index(pos);
        level{cnt}.X_val = zeros(size(pos));
        if cnt>1
            if ~is_cos
                level{cnt}.lowpass = fdct_wrapping_window((wedge_end(cnt)-level{cnt}.polar_r)/wedge_length(cnt));
            else
                level{cnt}.lowpass = gdct_cos((wedge_end(cnt)-level{cnt}.polar_r)/wedge_length(cnt));
            end
            level{cnt}.hipass = sqrt(1-(level{cnt}.lowpass).^2);
        end
        %temp(pos)= 10^10;
    end
    
    tiling = cell(1,nbscales-1);
    for cnt = 3:nbscales
        tiling{cnt-1} = cell(2,2);
        for cnt2 = 1:2
            if cnt2 == 1
                pos = find(level{cnt}.index>=0);
            else
                pos = find(level{cnt}.index<0);
            end
            tiling{cnt-1}{1,cnt2} = struct('index',0,'X_val',0);
            tiling{cnt-1}{1,cnt2}.index = level{cnt}.index(pos);
            tiling{cnt-1}{1,cnt2}.X_val = zeros(size(pos));
            if cnt2 == 1
                pos = find(level{cnt-1}.index>=0);
            else
                pos = find(level{cnt-1}.index<0);
            end
            tiling{cnt-1}{2,cnt2} = struct('index',0,'X_val',0);
            tiling{cnt-1}{2,cnt2}.index = level{cnt-1}.index(pos);
            tiling{cnt-1}{2,cnt2}.X_val = zeros(size(pos));
        end
    end
    tiling{1} = cell(1,1);
    
    %generate C
    fac = 1;
    temp = fftshift(fft(ifftshift(C{1}{1}{1})))/(sqrt(NG));
    x1 = mod(level{1}.index,NG)+1;
    XX(level{1}.index+pN+1) = temp(x1);
    x1 = mod(level{2}.index,NG)+1;
    XX(level{2}.index+pN+1) = temp(x1).*(level{2}.lowpass);
    if 1
        for cnt =3:nbscales
            if is_fac == 1
                fac = wedge_end(cnt)-wedge_end(cnt-1);
            else
                fac = 1;
            end
            for cnt2 = 1:2
                temp = fftshift(fft(ifftshift(C{1}{cnt-1}{cnt2})))/sqrt(NG)/fac;
                x1 = mod(tiling{cnt-1}{1,cnt2}.index+pN+1,NG)+1;
                if cnt2 == 1
                    pos = find(level{cnt}.index>=0);
                else
                    pos = find(level{cnt}.index<0);
                end
                tiling{cnt-1}{1,cnt2}.X_val = tiling{cnt-1}{1,cnt2}.X_val + temp(x1).*(level{cnt}.lowpass(pos));
                x1 = mod(tiling{cnt-1}{2,cnt2}.index+pN+1,NG)+1;
                if cnt2 == 1
                    pos = find(level{cnt-1}.index>=0);
                else
                    pos = find(level{cnt-1}.index<0);
                end
                tiling{cnt-1}{2,cnt2}.X_val = tiling{cnt-1}{2,cnt2}.X_val + temp(x1).*(level{cnt-1}.hipass(pos));
            end
        end
    end
    for cnt = 3:nbscales
        for cnt2 = 1:2
            if cnt2 == 1
                pos = find(level{cnt}.index>=0);
            else
                pos = find(level{cnt}.index<0);
            end
            level{cnt}.X_val(pos) = level{cnt}.X_val(pos) + tiling{cnt-1}{1,cnt2}.X_val;
            if cnt2 == 1
                pos = find(level{cnt-1}.index>=0);
            else
                pos = find(level{cnt-1}.index<0);
            end
            level{cnt-1}.X_val(pos) = level{cnt-1}.X_val(pos) + tiling{cnt-1}{2,cnt2}.X_val;
        end
    end
    for cnt = 2:nbscales
        XX(level{cnt}.index+pN+1) = XX(level{cnt}.index+pN+1) + level{cnt}.X_val;
    end
    %XX = XX(end:-1:1);%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if is_real
        if mod(N,2) == 0
            vec1 = 2:N;
        else
            vec1 = 1:N;
        end
        temp = XX(vec1);
        XX(vec1) = real(temp + temp(end:-1:1)) + i*imag(temp - temp(end:-1:1));
        
        if mod(N,2) == 0
            pos1 = N/2+1;
        else
            pos1 = (N+1)/2;
        end
        XX(pos1) = XX(pos1)/2;
    end
    x = fftshift(ifft2(ifftshift(XX)))*sqrt(N);
    amplt = abs(x);
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
    
    temp = wedge_end(cut_num);
    wedge_length = wedge_length(cut_num+2:end);
    wedge_length = [R_low-min(2,R_low-temp), min(2,R_low-temp), wedge_end(cut_num+1)-R_low, wedge_length];
    wedge_end = wedge_end(cut_num+1:end);
    wedge_end = [R_low-min(2,R_low-temp), R_low, wedge_end];
    
    nbscales = length(wedge_end);
    temp = cell(1,nbscales-1);
    for j = 2:nbscales-1
        temp{j} = C{1}{j-1};
    end
    C{1} = temp;
    
    %Polar coordinate
    Pos_radius = abs(Pos_index);
    
    level = cell(1,nbscales);%the information of fictitious tiling is nonnecessary
    %temp = Pos_radius;
    for cnt = 1:nbscales
        level{cnt} = struct('polar_r',0,'index',0,'X_val',0,'lowpass',0,'hipass',0);
        if cnt == 1
            temp = floor(wedge_end(cnt))+1;
            st = (-temp+1+(N-mod(N,2))/2);
            ed = (temp-2+(N+2-mod(N,2))/2);
            pos = st:ed;%find(temp<=wedge_end(cnt));
        else %if cnt<nbscales-1
            stOld = st;
            edOld = ed;
            temp = floor(wedge_end(cnt))+1;
            %st = max((-temp+2+(N-mod(N,2))/2),1);
            %ed = min((temp-1+(N+2-mod(N,2))/2),N);
            st = max((-temp+1+(N-mod(N,2))/2),1);
            ed = min((temp-2+(N+2-mod(N,2))/2),N);
            pos = [st:(stOld-1) (edOld+1):ed];%find(temp<=wedge_end(cnt));
        end
        level{cnt}.polar_r = Pos_radius(pos);
        level{cnt}.index = Pos_index(pos);
        level{cnt}.X_val = zeros(1,length(pos));
        if cnt>1
            if ~is_cos
                level{cnt}.lowpass = fdct_wrapping_window((wedge_end(cnt)-level{cnt}.polar_r)/wedge_length(cnt));
            else
                level{cnt}.lowpass = gdct_cos((wedge_end(cnt)-level{cnt}.polar_r)/wedge_length(cnt));
            end
            level{cnt}.hipass = sqrt(1-(level{cnt}.lowpass).^2);
        end
        %temp(pos)= 10^10;
    end
    
    tiling = cell(1,nbscales-1);
    for cnt = 3:nbscales
        tiling{cnt-1} = cell(2,2);
        for cnt2 = 1:2
            if cnt2 == 1
                pos = find(level{cnt}.index>=0);
            else
                pos = find(level{cnt}.index<0);
            end
            tiling{cnt-1}{1,cnt2} = struct('index',0,'X_val',0);
            tiling{cnt-1}{1,cnt2}.index = level{cnt}.index(pos);
            tiling{cnt-1}{1,cnt2}.X_val = zeros(1,length(pos));
            if cnt2 == 1
                pos = find(level{cnt-1}.index>=0);
            else
                pos = find(level{cnt-1}.index<0);
            end
            tiling{cnt-1}{2,cnt2} = struct('index',0,'X_val',0);
            tiling{cnt-1}{2,cnt2}.index = level{cnt-1}.index(pos);
            tiling{cnt-1}{2,cnt2}.X_val = zeros(1,length(pos));
        end
    end
    tiling{1} = cell(1,1);
    
    if 1
        for cnt =3:nbscales
            if is_fac == 1
                fac = wedge_end(cnt)-wedge_end(cnt-1);
            else
                fac = 1;
            end
            for cnt2 = 1:2
                temp = fftshift(fft(ifftshift(C{1}{cnt-1}{cnt2})))/sqrt(NG)/fac;
                x1 = mod(tiling{cnt-1}{1,cnt2}.index+pN+1,NG)+1;
                if cnt2 == 1
                    pos = find(level{cnt}.index>=0);
                else
                    pos = find(level{cnt}.index<0);
                end
                tiling{cnt-1}{1,cnt2}.X_val = tiling{cnt-1}{1,cnt2}.X_val + temp(x1).*(level{cnt}.lowpass(pos));
                x1 = mod(tiling{cnt-1}{2,cnt2}.index+pN+1,NG)+1;
                if cnt2 == 1
                    pos = find(level{cnt-1}.index>=0);
                else
                    pos = find(level{cnt-1}.index<0);
                end
                tiling{cnt-1}{2,cnt2}.X_val = tiling{cnt-1}{2,cnt2}.X_val + temp(x1).*(level{cnt-1}.hipass(pos));
            end
        end
    end
    for cnt = 3:nbscales
        for cnt2 = 1:2
            if cnt2 == 1
                pos = find(level{cnt}.index>=0);
            else
                pos = find(level{cnt}.index<0);
            end
            level{cnt}.X_val(pos) = level{cnt}.X_val(pos) + tiling{cnt-1}{1,cnt2}.X_val;
            if cnt2 == 1
                pos = find(level{cnt-1}.index>=0);
            else
                pos = find(level{cnt-1}.index<0);
            end
            level{cnt-1}.X_val(pos) = level{cnt-1}.X_val(pos) + tiling{cnt-1}{2,cnt2}.X_val;
        end
    end
    for cnt = 2:nbscales
        XX(level{cnt}.index+pN+1) = XX(level{cnt}.index+pN+1) + level{cnt}.X_val;
    end
    
    %XX = XX(end:-1:1,:);%!!!!!!!!!!!!!!!!!!!!!!!!!check
    if is_real
        if mod(N,2) == 0
            vec1 = 2:N;
        else
            vec1 = 1:N;
        end
        temp = XX(vec1);
        XX(vec1) = real(temp + temp(end:-1:1)) + i*imag(temp - temp(end:-1:1));
        
        if mod(N,2) == 0
            pos1 = N/2+1;
        else
            pos1 = (N+1)/2;
        end
        XX(pos1) = XX(pos1)/2;
    end
    x = fftshift(ifft2(ifftshift(XX)))*sqrt(N);
    amplt = abs(x);
    if is_real, x = real(x); end;
    
end
