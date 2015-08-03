function C = wp1_ext(x, is_real, is_unif, typeNUFFT, xo ,NG, R_high, R_low,  rad, is_cos, t_sc, red, is_fac, wedge_length_coarse)
% wp1_ext.m - derivative of 1D Wave Packet Transform in space variable
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
%   t_sc    scaling parameter for radius
%               [default set to 1-1/4]
%   red         redundancy parameter, red is a positive integer
%                [default set to 1]
%   is_fac      0: do not increase the magnitude of high frequency wave
%               packet coefficients; 1: increase;
%               [default set to 1, better to visualize high frequency
%               instantaneous frequencies]
%   wedge_length_coarse
%               length of coarsest wedge
%               [default set to 4]
%
% Outputs
%   C           Cell array of derivative of wave packet coefficients.
%               C{i}{j}{l}(k) is the coefficient at
%                   - the i'th frame, i=1,...,red
%                   - scale j: integer, from coarsest to finest scale,
%                   - l = 1 or 2,
%                   - position k: integers, size doesn't change
%
%by Haizhao Yang

if nargin < 3, is_unif = 1; end;
N = length(x);
if is_unif
    X = fftshift(fft(ifftshift(x)))/sqrt(N);
else
    error('Please use a non-uniform FFT!');
end

if mod(N,2) == 0
    pN = N/2;
    pos_vec = -pN:(pN-1);
else
    pN = (N-1)/2;
    pos_vec = -pN:pN;
end
Pos_index = pos_vec;

X = (2*pi*i*Pos_index).*X;

if nargin < 2, is_real = 0; end;
if is_real
    X( mod(N+1,2)+1: (N-mod(N,2))/2) = 0;
% else
%     R_low = 0;
end
if nargin < 10, is_cos = 1; end;
if nargin < 11, t_sc = 1 - 1/4; end;
if nargin < 12, red = 1; end;
if nargin < 13, is_fac = 1; end;
if nargin < 14, wedge_length_coarse = 4; end;
if wedge_length_coarse < 4
    wedge_length_coarse = 4;
    display('wedge_length_coarse is too small. Run with wedge_length_coarse = 4.');
end

% Initialization: data structure, decide wedge_length

R = (N-mod(N,2))/2;
if nargin > 6, R = min(R_high,R); end;
if nargin < 8, R_low = 0; end;
if ~is_real, R_low = 0; end;
if R_low <= wedge_length_coarse
    R_low = 0;
    display('R_low <= wedge_length_coarse, set R_low = 0.');
end
if nargin < 9, rad = 1; end;
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
NN = 0;
Pos_radius = abs(Pos_index);
if 0
    %Discard zero entries, can speed up
    pos = find(abs(X)>0);%can set up a threshold here !!!
    Pos_index = Pos_index(pos);
    Pos_radius = Pos_radius(pos);
    XX = X(pos);
else
    XX = X;
end
if R_low == 0
    nbscales = cell(1,red);
    C = cell(1,red);
    level = cell(1,red);
    tiling = cell(1,red);
    for cntred = 1:red
        nbscales{1,cntred} = length(wedge_end{1,cntred});
        C{1,cntred} = cell(1,nbscales{1,cntred}-1);
        for j = 1:nbscales{1,cntred}-1
            C{1,cntred}{j} = cell(2,1);
        end
        %check the sampling grid in the space domain
        NN = max(NN,ceil(wedge_length{1,cntred}(end-1)+wedge_length{1,cntred}(end-2)));
        level{cntred} = cell(1,nbscales{cntred});%the information of fictitious tiling is nonnecessary
        tiling{cntred} = cell(1,nbscales{cntred}-1);
    end
    if nargin < 7, NG = NN; display('Sampling grid not assigned. Run with the smallest grid.'); end
    if NG < NN
        display('Sampling grid in the space is too small! Run with the smallest grid.');
        NG = NN;
    end
    
    %temp2 = Pos_radius;
    for cntred = 1:red
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
            % maybe out of range
            level{cntred}{cnt}.polar_r = Pos_radius(pos);
            level{cntred}{cnt}.index = Pos_index(pos);
            level{cntred}{cnt}.X_val = XX(pos);
            if cnt>1
                if ~is_cos
                    level{cntred}{cnt}.lowpass = fdct_wrapping_window((wedge_end{cntred}(cnt)-level{cntred}{cnt}.polar_r)/wedge_length{cntred}(cnt));
                else
                    level{cntred}{cnt}.lowpass = gdct_cos((wedge_end{cntred}(cnt)-level{cntred}{cnt}.polar_r)/wedge_length{cntred}(cnt));
                end
                level{cntred}{cnt}.hipass = sqrt(1-(level{cntred}{cnt}.lowpass).^2);
            end
        end
        for cnt = 3:nbscales{cntred}
            tiling{cntred}{cnt-1} = cell(2,2);
            for cnt2 = 1:2
                % maybe out of range
                if cnt2 == 1
                    pos = find(level{cntred}{cnt}.index>=0);
                else
                    pos = find(level{cntred}{cnt}.index<0);
                end
                tiling{cntred}{cnt-1}{1,cnt2} = struct('index',0,'X_val',0);
                tiling{cntred}{cnt-1}{1,cnt2}.index = level{cntred}{cnt}.index(pos);
                tiling{cntred}{cnt-1}{1,cnt2}.X_val = level{cntred}{cnt}.X_val(pos).*(level{cntred}{cnt}.lowpass(pos));
                if cnt2 == 1
                    pos = find(level{cntred}{cnt-1}.index>=0);
                else
                    pos = find(level{cntred}{cnt-1}.index<0);
                end
                tiling{cntred}{cnt-1}{2,cnt2} = struct('index',0,'X_val',0);
                tiling{cntred}{cnt-1}{2,cnt2}.index = level{cntred}{cnt-1}.index(pos);
                tiling{cntred}{cnt-1}{2,cnt2}.X_val = level{cntred}{cnt-1}.X_val(pos).*(level{cntred}{cnt-1}.hipass(pos));
            end
        end
        tiling{cntred}{1} = cell(1,1);
        %generate C
        temp = zeros(1,NG);
        x1 = mod(level{cntred}{1}.index,NG)+1;
        temp(x1) = level{cntred}{1}.X_val;
        x1 = mod(level{cntred}{2}.index,NG)+1;
        temp(x1) = (level{cntred}{2}.X_val).*(level{cntred}{2}.lowpass);
        fac = 1;
        C{cntred}{1}{1} = fftshift(ifft(ifftshift(temp)))*sqrt(NG);
        for cnt =3:nbscales{cntred}
            if is_fac == 1
                fac = wedge_end{cntred}(cnt)-wedge_end{cntred}(cnt-1);
            else
                fac = 1;
            end
            for cnt2 = 1:2
                temp = zeros(1,NG);
                x1 = mod(tiling{cntred}{cnt-1}{1,cnt2}.index+pN+1,NG)+1;
                temp(x1) = tiling{cntred}{cnt-1}{1,cnt2}.X_val;
                x1 = mod(tiling{cntred}{cnt-1}{2,cnt2}.index+pN+1,NG)+1;
                temp(x1) = tiling{cntred}{cnt-1}{2,cnt2}.X_val;
                C{cntred}{cnt-1}{cnt2} = fftshift(ifft(ifftshift(temp)))*sqrt(NG)*fac;
            end
        end
    end
else
    %--------------------------------------------------------------------------
    %cut tiling
    %cut the tiling
    for cntred = 1:red
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
        
        C{cntred} = cell(1,nbscales{cntred}-1);
        for j = 1:nbscales{cntred}-1
            C{cntred}{j} = cell(2,1);
        end
        %check the sampling grid in the space domain
        NN = max(NN,ceil(wedge_length{cntred}(end-1)+wedge_length{cntred}(end-2)));
        level{cntred} = cell(1,nbscales{cntred});%the information of fictitious tiling is nonnecessary
        tiling{cntred} = cell(1,nbscales{cntred}-1);
    end
    if nargin < 7, NG = NN; display('Sampling grid not assigned. Run with the smallest grid.'); end
    if NG<NN
        display('Sampling grid in the space is too small! Run with the smallest grid.');
        NG = NN;
    end
    
    %temp = Pos_radius;
    for cntred = 1:red
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
            level{cntred}{cnt}.X_val = XX(pos);
            if cnt>1
                if ~is_cos
                    level{cntred}{cnt}.lowpass = fdct_wrapping_window((wedge_end{cntred}(cnt)-level{cntred}{cnt}.polar_r)/wedge_length{cntred}(cnt));
                else
                    level{cntred}{cnt}.lowpass = gdct_cos((wedge_end{cntred}(cnt)-level{cntred}{cnt}.polar_r)/wedge_length{cntred}(cnt));
                end
                level{cntred}{cnt}.hipass = sqrt(1-(level{cntred}{cnt}.lowpass).^2);
            end
            %temp(pos)= 10^10;
        end
        
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
                tiling{cntred}{cnt-1}{1,cnt2}.X_val = level{cntred}{cnt}.X_val(pos).*(level{cntred}{cnt}.lowpass(pos));
                if cnt2 == 1
                    pos = find(level{cntred}{cnt-1}.index>=0);
                else
                    pos = find(level{cntred}{cnt-1}.index<0);
                end
                tiling{cntred}{cnt-1}{2,cnt2} = struct('index',0,'X_val',0);
                tiling{cntred}{cnt-1}{2,cnt2}.index = level{cntred}{cnt-1}.index(pos);
                tiling{cntred}{cnt-1}{2,cnt2}.X_val = level{cntred}{cnt-1}.X_val(pos).*(level{cntred}{cnt-1}.hipass(pos));
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
                temp = zeros(1,NG);
                x1 = mod(tiling{cntred}{cnt-1}{1,cnt2}.index+pN+1,NG)+1;
                temp(x1) = tiling{cntred}{cnt-1}{1,cnt2}.X_val;
                x1 = mod(tiling{cntred}{cnt-1}{2,cnt2}.index+pN+1,NG)+1;
                temp(x1) = tiling{cntred}{cnt-1}{2,cnt2}.X_val;
                C{cntred}{cnt-1}{cnt2} = fftshift(ifft(ifftshift(temp)))*sqrt(NG)*fac;
            end
        end
        temp = cell(1,nbscales{cntred}-2);
        for cnt = 1:nbscales{cntred}-2
            temp{cnt} = C{cntred}{cnt+1};
        end
        C{cntred} = temp;
    end
end