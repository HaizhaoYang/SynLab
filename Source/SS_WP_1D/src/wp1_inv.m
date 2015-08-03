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
red = length(C);
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
wedge_length = cell(1,red);
wedge_end = cell(1,red);
da = rad*R^t_sc/2/red;
for cntred = 1:red
    wedge_length{1,cntred} = 2;
    wedge_end{1,cntred} = [R R+2]+(cntred-1)*da;
    while wedge_end{1,cntred}(1) > wedge_length_coarse | length(wedge_end{1,cntred})<4
        Length = rad*wedge_end{1,cntred}(1)^t_sc/2;
        wedge_length{1,cntred} = [Length wedge_length{1,cntred}];
        wedge_end{1,cntred} = [(wedge_end{1,cntred}(1)-Length) wedge_end{1,cntred}];
    end
    wedge_length{1,cntred} = [wedge_end{1,cntred}(1) wedge_length{1,cntred}];
end

if mod(N,2) == 0
    pN = N/2;
    pos_vec = -pN:(pN-1);
else
    pN = (N-1)/2;
    pos_vec = -pN:pN;
end
Pos_index = pos_vec;

if R_low == 0
    XX = cell(1,red);
    nbscales = cell(1,red);
    level = cell(1,red);
    Pos_radius = abs(Pos_index);
    for cntred = 1:red
        XX{cntred} = zeros(1,N);
        nbscales{cntred} = length(wedge_end{cntred});
        level{cntred} = cell(1,nbscales{cntred});%the information of fictitious tiling is nonnecessary
        for cnt = 1:nbscales{cntred}
            level{cntred}{cnt} = struct('polar_r',0,'index',0,'X_val',0,'lowpass',0,'hipass',0);
            if cnt == 1
                temp = floor(wedge_end{cntred}(cnt))+1;
                %st = (-temp+2+(N-mod(N,2))/2);
                %ed = (temp-1+(N+2-mod(N,2))/2);
                st = (-temp+1+(N-mod(N,2))/2);
                ed = (temp-2+(N+2-mod(N,2))/2);
                pos = st:ed;%find(temp<=wedge_end(cnt));
            else %if cnt<nbscales-1
                stOld = st;
                edOld = ed;
                temp = floor(wedge_end{cntred}(cnt))+1;
                %st = max((-temp+2+(N-mod(N,2))/2),1);
                %ed = min((temp-1+(N+2-mod(N,2))/2),N);
                st = max((-temp+1+(N-mod(N,2))/2),1);
                ed = min((temp-2+(N+2-mod(N,2))/2),N);
                pos = [st:(stOld-1) (edOld+1):ed];%find(temp<=wedge_end(cnt));
            end
            level{cntred}{cnt}.polar_r = Pos_radius(pos);
            level{cntred}{cnt}.index = Pos_index(pos);
            level{cntred}{cnt}.X_val = zeros(size(pos));
            if cnt>1
                if ~is_cos
                    level{cntred}{cnt}.lowpass = fdct_wrapping_window((wedge_end{cntred}(cnt)-level{cntred}{cnt}.polar_r)/wedge_length{cntred}(cnt));
                else
                    level{cntred}{cnt}.lowpass = gdct_cos((wedge_end{cntred}(cnt)-level{cntred}{cnt}.polar_r)/wedge_length{cntred}(cnt));
                end
                level{cntred}{cnt}.hipass = sqrt(1-(level{cntred}{cnt}.lowpass).^2);
            end
        end
    end
    tiling = cell(1,red);
    for cntred = 1:red
        tiling{cntred} = cell(1,nbscales{cntred}-1);
        for cnt = 3:nbscales{cntred}
            tiling{cntred}{cnt-1} = cell(2,2);
            for cnt2 = 1:2
                if cnt2 == 1
                    pos = find(level{cntred}{cnt}.index>=0);
                else
                    pos = find(level{cntred}{cnt}.index<0);
                end
                tiling{cntred}{cnt-1}{1,cnt2} = struct('index',0,'X_val',0);
                tiling{cntred}{cnt-1}{1,cnt2}.index = level{cntred}{cnt}.index(pos);
                tiling{cntred}{cnt-1}{1,cnt2}.X_val = zeros(size(pos));
                if cnt2 == 1
                    pos = find(level{cntred}{cnt-1}.index>=0);
                else
                    pos = find(level{cntred}{cnt-1}.index<0);
                end
                tiling{cntred}{cnt-1}{2,cnt2} = struct('index',0,'X_val',0);
                tiling{cntred}{cnt-1}{2,cnt2}.index = level{cntred}{cnt-1}.index(pos);
                tiling{cntred}{cnt-1}{2,cnt2}.X_val = zeros(size(pos));
            end
        end
        tiling{cntred}{1} = cell(1,1);
    end
    for cntred = 1:red
        temp = fftshift(fft(ifftshift(C{cntred}{1}{1})))/(sqrt(NG));
        x1 = mod(level{cntred}{1}.index,NG)+1;
        XX{cntred}(level{cntred}{1}.index+pN+1) = temp(x1);
        x1 = mod(level{cntred}{2}.index,NG)+1;
        XX{cntred}(level{cntred}{2}.index+pN+1) = temp(x1).*(level{cntred}{2}.lowpass);
        for cnt =3:nbscales{cntred}
            if is_fac == 1
                fac = wedge_end{cntred}(cnt)-wedge_end{cntred}(cnt-1);
            else
                fac = 1;
            end
            for cnt2 = 1:2
                temp = fftshift(fft(ifftshift(C{cntred}{cnt-1}{cnt2})))/sqrt(NG)/fac;
                x1 = mod(tiling{cntred}{cnt-1}{1,cnt2}.index+pN+1,NG)+1;
                if cnt2 == 1
                    pos = find(level{cntred}{cnt}.index>=0);
                else
                    pos = find(level{cntred}{cnt}.index<0);
                end
                tiling{cntred}{cnt-1}{1,cnt2}.X_val = tiling{cntred}{cnt-1}{1,cnt2}.X_val + temp(x1).*(level{cntred}{cnt}.lowpass(pos));
                x1 = mod(tiling{cntred}{cnt-1}{2,cnt2}.index+pN+1,NG)+1;
                if cnt2 == 1
                    pos = find(level{cntred}{cnt-1}.index>=0);
                else
                    pos = find(level{cntred}{cnt-1}.index<0);
                end
                tiling{cntred}{cnt-1}{2,cnt2}.X_val = tiling{cntred}{cnt-1}{2,cnt2}.X_val + temp(x1).*(level{cntred}{cnt-1}.hipass(pos));
            end
        end
        for cnt = 3:nbscales{cntred}
            for cnt2 = 1:2
                if cnt2 == 1
                    pos = find(level{cntred}{cnt}.index>=0);
                else
                    pos = find(level{cntred}{cnt}.index<0);
                end
                level{cntred}{cnt}.X_val(pos) = level{cntred}{cnt}.X_val(pos) + tiling{cntred}{cnt-1}{1,cnt2}.X_val;
                if cnt2 == 1
                    pos = find(level{cntred}{cnt-1}.index>=0);
                else
                    pos = find(level{cntred}{cnt-1}.index<0);
                end
                level{cntred}{cnt-1}.X_val(pos) = level{cntred}{cnt-1}.X_val(pos) + tiling{cntred}{cnt-1}{2,cnt2}.X_val;
            end
        end
        for cnt = 2:nbscales{cntred}
            XX{cntred}(level{cntred}{cnt}.index+pN+1) = XX{cntred}(level{cntred}{cnt}.index+pN+1) + level{cntred}{cnt}.X_val;
        end
        %XX = XX(end:-1:1);%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if is_real
            if mod(N,2) == 0
                vec1 = 2:N;
            else
                vec1 = 1:N;
            end
            temp = XX{cntred}(vec1);
            XX{cntred}(vec1) = real(temp + temp(end:-1:1)) + i*imag(temp - temp(end:-1:1));
            
            if mod(N,2) == 0
                pos1 = N/2+1;
            else
                pos1 = (N+1)/2;
            end
            XX{cntred}(pos1) = XX{cntred}(pos1)/2;
        end
    end
    x = zeros(1,N);
    for cntred = 1:red
        x = x + fftshift(ifft2(ifftshift(XX{cntred})))*sqrt(N);
    end
    x = x/red;
    amplt = abs(x);
    if is_real, x = real(x); end;
    
else
    XX = cell(1,red);
    Pos_radius = abs(Pos_index);
    level = cell(1,red);
    tiling = cell(1,red);
    for cntred = 1:red
        XX{cntred} = zeros(1,N);
        cut_pos = find(wedge_end{cntred}<R_low);
        cut_num = length(cut_pos);
        while cut_num < 2
            Length = rad*wedge_end{cntred}(1)^t_sc/2;
            wedge_end{cntred} = [(wedge_end{cntred}(1)-Length) wedge_end{cntred}];
            wedge_length{cntred} = wedge_length{cntred}(2:end);
            wedge_length{cntred} = [wedge_end{cntred}(1) Length wedge_length{cntred}];
            cut_pos = find(wedge_end{cntred}<R_low);
            cut_num = length(cut_pos);
        end
        
        
        temp = wedge_end{cntred}(cut_num);
        wedge_length{cntred} = wedge_length{cntred}(cut_num+2:end);
        wedge_length{cntred} = [R_low-min(2,R_low-temp), min(2,R_low-temp), wedge_end{cntred}(cut_num+1)-R_low, wedge_length{cntred}];
        wedge_end{cntred} = wedge_end{cntred}(cut_num+1:end);
        wedge_end{cntred} = [R_low-min(2,R_low-temp), R_low, wedge_end{cntred}];
        
        nbscales{cntred} = length(wedge_end{cntred});
        temp = cell(1,nbscales{cntred}-1);
        for j = 2:nbscales{cntred}-1
            temp{j} = C{cntred}{j-1};
        end
        C{cntred} = temp;
    end
    
    for cntred = 1:red
        level{cntred} = cell(1,nbscales{cntred});%the information of fictitious tiling is nonnecessary
        for cnt = 1:nbscales{cntred}
            level{cntred}{cnt} = struct('polar_r',0,'index',0,'X_val',0,'lowpass',0,'hipass',0);
            if cnt == 1
                temp = floor(wedge_end{cntred}(cnt))+1;
                st = (-temp+1+(N-mod(N,2))/2);
                ed = (temp-2+(N+2-mod(N,2))/2);
                pos = st:ed;%find(temp<=wedge_end(cnt));
            else %if cnt<nbscales-1
                stOld = st;
                edOld = ed;
                temp = floor(wedge_end{cntred}(cnt))+1;
                %st = max((-temp+2+(N-mod(N,2))/2),1);
                %ed = min((temp-1+(N+2-mod(N,2))/2),N);
                st = max((-temp+1+(N-mod(N,2))/2),1);
                ed = min((temp-2+(N+2-mod(N,2))/2),N);
                pos = [st:(stOld-1) (edOld+1):ed];%find(temp<=wedge_end(cnt));
            end
            level{cntred}{cnt}.polar_r = Pos_radius(pos);
            level{cntred}{cnt}.index = Pos_index(pos);
            level{cntred}{cnt}.X_val = zeros(1,length(pos));
            if cnt>1
                if ~is_cos
                    level{cntred}{cnt}.lowpass = fdct_wrapping_window((wedge_end{cntred}(cnt)-level{cntred}{cnt}.polar_r)/wedge_length{cntred}(cnt));
                else
                    level{cntred}{cnt}.lowpass = gdct_cos((wedge_end{cntred}(cnt)-level{cntred}{cnt}.polar_r)/wedge_length{cntred}(cnt));
                end
                level{cntred}{cnt}.hipass = sqrt(1-(level{cntred}{cnt}.lowpass).^2);
            end
        end
    end
    
    for cntred = 1:red
        tiling{cntred} = cell(1,nbscales{cntred}-1);
        for cnt = 3:nbscales{cntred}
            tiling{cntred}{cnt-1} = cell(2,2);
            for cnt2 = 1:2
                if cnt2 == 1
                    pos = find(level{cntred}{cnt}.index>=0);
                else
                    pos = find(level{cntred}{cnt}.index<0);
                end
                tiling{cntred}{cnt-1}{1,cnt2} = struct('index',0,'X_val',0);
                tiling{cntred}{cnt-1}{1,cnt2}.index = level{cntred}{cnt}.index(pos);
                tiling{cntred}{cnt-1}{1,cnt2}.X_val = zeros(1,length(pos));
                if cnt2 == 1
                    pos = find(level{cntred}{cnt-1}.index>=0);
                else
                    pos = find(level{cntred}{cnt-1}.index<0);
                end
                tiling{cntred}{cnt-1}{2,cnt2} = struct('index',0,'X_val',0);
                tiling{cntred}{cnt-1}{2,cnt2}.index = level{cntred}{cnt-1}.index(pos);
                tiling{cntred}{cnt-1}{2,cnt2}.X_val = zeros(1,length(pos));
            end
        end
        tiling{cntred}{1} = cell(1,1);
        
        for cnt =3:nbscales{cntred}
            if is_fac == 1
                fac = wedge_end{cntred}(cnt)-wedge_end{cntred}(cnt-1);
            else
                fac = 1;
            end
            for cnt2 = 1:2
                temp = fftshift(fft(ifftshift(C{cntred}{cnt-1}{cnt2})))/sqrt(NG)/fac;
                x1 = mod(tiling{cntred}{cnt-1}{1,cnt2}.index+pN+1,NG)+1;
                if cnt2 == 1
                    pos = find(level{cntred}{cnt}.index>=0);
                else
                    pos = find(level{cntred}{cnt}.index<0);
                end
                if length(pos)>0
                    tiling{cntred}{cnt-1}{1,cnt2}.X_val = tiling{cntred}{cnt-1}{1,cnt2}.X_val + temp(x1).*(level{cntred}{cnt}.lowpass(pos));
                end
                x1 = mod(tiling{cntred}{cnt-1}{2,cnt2}.index+pN+1,NG)+1;
                if cnt2 == 1
                    pos = find(level{cntred}{cnt-1}.index>=0);
                else
                    pos = find(level{cntred}{cnt-1}.index<0);
                end
                if length(pos)>0
                    tiling{cntred}{cnt-1}{2,cnt2}.X_val = tiling{cntred}{cnt-1}{2,cnt2}.X_val + temp(x1).*(level{cntred}{cnt-1}.hipass(pos));
                end
            end
        end
        for cnt = 3:nbscales{cntred}
            for cnt2 = 1:2
                if cnt2 == 1
                    pos = find(level{cntred}{cnt}.index>=0);
                else
                    pos = find(level{cntred}{cnt}.index<0);
                end
                if length(pos)>0
                    level{cntred}{cnt}.X_val(pos) = level{cntred}{cnt}.X_val(pos) + tiling{cntred}{cnt-1}{1,cnt2}.X_val;
                end
                if cnt2 == 1
                    pos = find(level{cntred}{cnt-1}.index>=0);
                else
                    pos = find(level{cntred}{cnt-1}.index<0);
                end
                if length(pos)>0
                    level{cntred}{cnt-1}.X_val(pos) = level{cntred}{cnt-1}.X_val(pos) + tiling{cntred}{cnt-1}{2,cnt2}.X_val;
                end
            end
        end
        for cnt = 2:nbscales{cntred}
            XX{cntred}(level{cntred}{cnt}.index+pN+1) = XX{cntred}(level{cntred}{cnt}.index+pN+1) + level{cntred}{cnt}.X_val;
        end
        
        %XX = XX(end:-1:1,:);%!!!!!!!!!!!!!!!!!!!!!!!!!check
        if is_real
            if mod(N,2) == 0
                vec1 = 2:N;
            else
                vec1 = 1:N;
            end
            temp = XX{cntred}(vec1);
            XX{cntred}(vec1) = real(temp + temp(end:-1:1)) + i*imag(temp - temp(end:-1:1));
            
            if mod(N,2) == 0
                pos1 = N/2+1;
            else
                pos1 = (N+1)/2;
            end
            XX{cntred}(pos1) = XX{cntred}(pos1)/2;
        end
    end
    x = zeros(1,N);
    for cntred = 1:red
        x = x + fftshift(ifft2(ifftshift(XX{cntred})))*sqrt(N);
    end
    x = x/red;
    amplt = abs(x);
    if is_real, x = real(x); end;
    
end
