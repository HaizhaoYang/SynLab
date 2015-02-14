function [kb avgdx avgdy] = SS_ct_v2(num_direction,fff,sz,NB,rad,is_real,radiusRange,angleRange,epsl,red,t_sc, s_sc,is_cos,  wedge_length_coarse)
% This code provides the synchrosqueezed energy in a given range of
% frequency. The energy is given in a Cartesian coordinate.
%
% Input:
% fff is the image
% sz(1) is the number of samples in the vertical direction in the image
% sz(2) is the number of samples in the horizontal direction in the image
% (2*NB+1) by (2*NB+1) is the size of samples in the phase in [-R_high,R_high]*[-R_high,R_high]
% rad is the smallest size of supports of wave packets
% is_real = 0: fff is complex
% is_real = 1: fff is real
% radiusRange = [R_low R_high], where
% R_low: lower bound of interested frequency
% R_high: upper bound of interested frequency
% angleRange = [angleRange(1) angleRange(2)], where
% angleRange(1): lower bound of interested angles
% angleRange(2): upper bound of interested angles
% epsl: the threshold for small coefficients of general curvelet transform
% red is the redundancy parameter, red(1) for angle, red(2) for radii
% default set to [1 1]
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
% Output:
% kb the synchrosqueezed energy distribution in a Cartesian coordinate
% avgdx the weighted position in x axis
% avgdy the weighted position in y axis
% avgdx and avgdy are used to compute the weighted local wave vector
% estimate
%
%by Haizhao Yang

%figure;hold on;
%tmpf = zeros(size(fff));
%% set up parameters
%==========================================================================
if nargin < 5, rad = 1; end;
if rad > 2, display('rad too large! rad is set up to be 2'); rad = 2; end;
if nargin < 8, angleRange = [0 2*pi]; end;
temp = floor(angleRange(1)/2/pi);
angleRange = angleRange - 2*pi*temp;
if angleRange(1)>=angleRange(2) | angleRange(2)-angleRange(1)>2*pi, error('Invalid angleRange'); end;
if nargin < 9 , epsl = 1e-4; end;
if nargin < 6, is_real = 0; end;
if nargin < 10, red = [1 1]; end;
if nargin < 11, t_sc = 1 - 1/8; end;
if nargin < 12, s_sc = 1/2 + 1/8; end;
if nargin < 13, is_cos = 1; end;
if nargin < 14, wedge_length_coarse = 8; end;
if wedge_length_coarse < 4
    wedge_length_coarse = 4;
    display('wedge_length_coarse is too small. Run with wedge_length_coarse = 4.');
end

%set up size of image
[N1 N2] = size(fff);
N = max(N1,N2);

%% apply the general curvelet transform and squeeze nonzero energy
%==========================================================================
X = fftshift(fft2(ifftshift(fff)))/sqrt(N1*N2);
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
X1 = (2*pi*i*Pos_index_1).*X;
X2 = (2*pi*i*Pos_index_2).*X;

if is_real
    X( mod(N1+1,2)+1:end , mod(N2+1,2)+1: (N2-mod(N2,2))/2) = 0;
    X( mod(N1+1,2)+1: (N1-mod(N1,2))/2, (N2+1+mod(N2+1,2))/2) = 0;
    X1( mod(N1+1,2)+1:end , mod(N2+1,2)+1: (N2-mod(N2,2))/2) = 0;
    X1( mod(N1+1,2)+1: (N1-mod(N1,2))/2, (N2+1+mod(N2+1,2))/2) = 0;
    X2( mod(N1+1,2)+1:end , mod(N2+1,2)+1: (N2-mod(N2,2))/2) = 0;
    X2( mod(N1+1,2)+1: (N1-mod(N1,2))/2, (N2+1+mod(N2+1,2))/2) = 0;
end

% Initialization: data structure, decide wedge_length and wedge_end
R = sqrt(N1^2+N2^2)/2;
if nargin > 7, R_low = radiusRange(1); R_high = radiusRange(2); R = min(R_high,R); end;
if nargin < 7, R_low = 0;  R_high = R; end;
if R_low <= wedge_length_coarse
    R_low = 0;
    display('R_low <= wedge_length_coarse, set R_low = 0.');
end
if R<=R_low, error('R_low is too large! R_low at most %f',R); end;

%Polar coordinate
Pos_radius = sqrt(Pos_index_1.^2+Pos_index_2.^2);
pos = find(Pos_radius>0);
Pos_angle = zeros(N1,N2);
Pos_angle(pos) = acos(Pos_index_1(pos)./Pos_radius(pos));
pos = find(Pos_index_2<0);
Pos_angle(pos) = 2*pi-Pos_angle(pos);
if 1
    %Discard zero entries, can speed up
    pos = find(abs(X)>0);%can set up a threshold here !!!
    Pos_index_1 = Pos_index_1(pos);
    Pos_index_2 = Pos_index_2(pos);
    Pos_radius = Pos_radius(pos);
    Pos_angle = Pos_angle(pos);
    XX = X(pos);
    XX1 = X1(pos);
    XX2 = X2(pos);
    clear X X1 X2;
else
    XX = X(:);
    XX1 = X1(:);
    XX2 = X2(:);
    clear X X1 X2;
end


wedge_length_coarse = wedge_length_coarse*rad;
wedge_length = cell(1,red(2));
nbangles = cell(1,red(2));
angles = cell(1,red(2));
nbscales = cell(1,red(2));
da = rad*R^t_sc/2/red(2);
NN = 0;
for cntred = 1:red(2)
    wedge_length{1,cntred} = 2;
    wedge_end{1,cntred} = [R R+2]+(cntred-1)*da;
    while wedge_end{1,cntred}(1) > wedge_length_coarse | length(wedge_end{1,cntred})<4
        Length = rad*wedge_end{1,cntred}(1)^t_sc/2;
        wedge_length{1,cntred} = [Length wedge_length{1,cntred}];
        wedge_end{1,cntred} = [(wedge_end{1,cntred}(1)-Length) wedge_end{1,cntred}];
    end
    wedge_length{1,cntred} = [wedge_end{1,cntred}(1) wedge_length{1,cntred}];
    if R_low>0
        cut_pos = find(wedge_end{1,cntred}<R_low);
        cut_num = length(cut_pos);
        while cut_num < 2
            Length = rad*wedge_end{1,cntred}(1)^t_sc/2;
            wedge_end{1,cntred} = [(wedge_end{1,cntred}(1)-Length) wedge_end{1,cntred}];
            wedge_length{1,cntred} = wedge_length{1,cntred}(2:end);
            wedge_length{1,cntred} = [wedge_end{1,cntred}(1) Length wedge_length{1,cntred}];
            cut_pos = find(wedge_end{1,cntred}<R_low);
            cut_num = length(cut_pos);
        end
        temp = wedge_length{1,cntred}(3:end) + wedge_length{1,cntred}(2:end-1);
        arc_length = temp.^(s_sc/t_sc);
        nbangles{1,cntred} = round(pi*wedge_end{1,cntred}(3:end)./arc_length/2)*4;
        nbangles{1,cntred} = [1 1 nbangles{1,cntred}];
        
        temp = wedge_end{1,cntred}(cut_num);
        wedge_length{1,cntred} = wedge_length{1,cntred}(cut_num+2:end);
        wedge_length{1,cntred} = [R_low-min(2,R_low-temp), min(2,R_low-temp), wedge_end{1,cntred}(cut_num+1)-R_low, wedge_length{1,cntred}];
        wedge_end{1,cntred} = wedge_end{1,cntred}(cut_num+1:end);
        wedge_end{1,cntred} = [R_low-min(2,R_low-temp), R_low, wedge_end{1,cntred}];
        
        nbangles{1,cntred} = nbangles{1,cntred}(cut_num+1:end);
        nbangles{1,cntred} = [1 1 nbangles{1,cntred}];
        angles{1,cntred} = 2*pi./nbangles{1,cntred};
    else
        temp = wedge_length{1,cntred}(3:end) + wedge_length{1,cntred}(2:end-1);
        arc_length = temp.^(s_sc/t_sc);
        nbangles{1,cntred} = round(pi*wedge_end{1,cntred}(3:end)./arc_length/2)*4;
        nbangles{1,cntred} = [1 1 nbangles{1,cntred}];
        angles{1,cntred} = 2*pi./nbangles{1,cntred};
    end
    nbscales{1,cntred} = length(wedge_end{1,cntred});
    %check the sampling grid in the space domain
    NN = max(NN,ceil(max([wedge_length{1,cntred}(end-1)+wedge_length{1,cntred}(end-2),sqrt(wedge_end{1,cntred}(end-1)^2+wedge_end{1,cntred}(end-3)^2-2*wedge_end{1,cntred}(end-3)*wedge_end{1,cntred}(end-1)*cos(angles{1,cntred}(end-1))), 2*wedge_end{1,cntred}(end-1)*sin(angles{1,cntred}(end-1)) ])));
end
NN = ceil(NN/4)*4;
if nargin < 3, sz = [NN NN]; display('Sampling grid not assigned. Run with the smallest grid.'); end
pos = find(sz<NN);
if pos
    %display('Sampling grid in the space is too small! Run with the smallest grid.');
    sz(pos) = NN;
end
EXT = 10^10;
%dr = (R_high-R_low)/(NB(1));
%da = pi/(NB(2))/num_direction;
% kb = zeros(NB(1)+1,NB(2)+1,sz(1),sz(2));
% avgdx = zeros(NB(1)+1,NB(2)+1,sz(1),sz(2));
% avgdy = zeros(NB(1)+1,NB(2)+1,sz(1),sz(2));

WB = NB/R_high;%2*min(N,2*R_high)/NB;
kb = zeros(2*NB+1,2*NB+1,sz(1),sz(2));
avgdx = zeros(2*NB+1,2*NB+1,sz(1),sz(2));
avgdy = zeros(2*NB+1,2*NB+1,sz(1),sz(2));

for cntredR = 1:red(2)
    level = cell(1,nbscales{1,cntredR});%the information of fictitious tiling is nonnecessary
    temp = Pos_radius;
    for cnt = 1:nbscales{1,cntredR}
        level{cnt} = struct('polar_r',0,'polar_a',0,'index_1',0,'index_2',0,'X_val',0,'X1_val',0,'X2_val',0,'lowpass',0,'hipass',0);
        pos = find(temp<=wedge_end{1,cntredR}(cnt));
        level{cnt}.polar_r = Pos_radius(pos);
        level{cnt}.polar_a = Pos_angle(pos);
        level{cnt}.index_1 = Pos_index_1(pos);
        level{cnt}.index_2 = Pos_index_2(pos);
        level{cnt}.X_val = XX(pos);
        level{cnt}.X1_val = XX1(pos);
        level{cnt}.X2_val = XX2(pos);
        if cnt>1
            if ~is_cos
                level{cnt}.lowpass = fdct_wrapping_window((wedge_end{1,cntredR}(cnt)-level{cnt}.polar_r)/wedge_length{1,cntredR}(cnt));
            else
                level{cnt}.lowpass = gdct_cos((wedge_end{1,cntredR}(cnt)-level{cnt}.polar_r)/wedge_length{1,cntredR}(cnt));
            end
            level{cnt}.hipass = sqrt(1-(level{cnt}.lowpass).^2);
        end
        temp(pos)= 10^10;
    end
    %clear XX XX1 XX2;
    tiling = cell(1,nbscales{1,cntredR}-1);
    for cntredA = 0:red(1)-1
        for cnt = 3:nbscales{1,cntredR}
            tiling{cnt-1} = cell(2,nbangles{1,cntredR}(cnt));
            if angleRange(1)+2*pi-angleRange(2)<2*angles{1,cntredR}(cnt)
                cntRange = 1:nbangles{1,cntredR}(cnt);
            else
                cnt_st = mod(floor(angleRange(1)/angles{1,cntredR}(cnt)-cntredA/red(1))-1-1,nbangles{1,cntredR}(cnt))+1;
                cnt_ed = mod(ceil(angleRange(2)/angles{1,cntredR}(cnt)-cntredA/red(1))-1+1,nbangles{1,cntredR}(cnt))+1;
                if cnt_st<=cnt_ed
                    cntRange = cnt_st:cnt_ed;
                else
                    cntRange = [1:cnt_ed cnt_st:nbangles{1,cntredR}(cnt)];
                end
            end
            for cnt2 = cntRange
                %for cnt2 = 1:nbangles{1,cntredR}(cnt)
                if cnt2 < nbangles{1,cntredR}(cnt)
                    pos = find(level{cnt}.polar_a>=(cnt2-1+cntredA/red(1))*angles{1,cntredR}(cnt) & level{cnt}.polar_a < (cnt2+cntredA/red(1))*angles{1,cntredR}(cnt));
                    tiling{cnt-1}{1,cnt2} = struct('polar_a',0,'index_1',0,'index_2',0,'X_val',0,'X1_val',0,'X2_val',0,'lowpass',0,'hipass',0);
                    tiling{cnt-1}{1,cnt2}.polar_a = level{cnt}.polar_a(pos);
                    tiling{cnt-1}{1,cnt2}.index_1 = level{cnt}.index_1(pos);
                    tiling{cnt-1}{1,cnt2}.index_2 = level{cnt}.index_2(pos);
                    tiling{cnt-1}{1,cnt2}.X_val = level{cnt}.X_val(pos);
                    tiling{cnt-1}{1,cnt2}.X1_val = level{cnt}.X1_val(pos);
                    tiling{cnt-1}{1,cnt2}.X2_val = level{cnt}.X2_val(pos);
                    if ~is_cos
                        tiling{cnt-1}{1,cnt2}.lowpass = fdct_wrapping_window((tiling{cnt-1}{1,cnt2}.polar_a-(cnt2-1+cntredA/red(1))*angles{1,cntredR}(cnt))/angles{1,cntredR}(cnt));
                    else
                        tiling{cnt-1}{1,cnt2}.lowpass = gdct_cos((tiling{cnt-1}{1,cnt2}.polar_a-(cnt2-1+cntredA/red(1))*angles{1,cntredR}(cnt))/angles{1,cntredR}(cnt));
                    end
                    tiling{cnt-1}{1,cnt2}.hipass = sqrt(1-(tiling{cnt-1}{1,cnt2}.lowpass).^2);
                    tiling{cnt-1}{1,cnt2}.lowpass = (tiling{cnt-1}{1,cnt2}.lowpass).*(level{cnt}.lowpass(pos));
                    tiling{cnt-1}{1,cnt2}.hipass = (tiling{cnt-1}{1,cnt2}.hipass).*(level{cnt}.lowpass(pos));
                    pos = find(level{cnt-1}.polar_a>=(cnt2-1+cntredA/red(1))*angles{1,cntredR}(cnt) & level{cnt-1}.polar_a<(cnt2+cntredA/red(1))*angles{1,cntredR}(cnt));
                    tiling{cnt-1}{2,cnt2} = struct('polar_a',0,'index_1',0,'index_2',0,'X_val',0,'X1_val',0,'X2_val',0,'lowpass',0,'hipass',0);
                    tiling{cnt-1}{2,cnt2}.polar_a = level{cnt-1}.polar_a(pos);
                    tiling{cnt-1}{2,cnt2}.index_1 = level{cnt-1}.index_1(pos);
                    tiling{cnt-1}{2,cnt2}.index_2 = level{cnt-1}.index_2(pos);
                    tiling{cnt-1}{2,cnt2}.X_val = level{cnt-1}.X_val(pos);
                    tiling{cnt-1}{2,cnt2}.X1_val = level{cnt-1}.X1_val(pos);
                    tiling{cnt-1}{2,cnt2}.X2_val = level{cnt-1}.X2_val(pos);
                    if ~is_cos
                        tiling{cnt-1}{2,cnt2}.lowpass = fdct_wrapping_window((tiling{cnt-1}{2,cnt2}.polar_a-(cnt2-1+cntredA/red(1))*angles{1,cntredR}(cnt))/angles{1,cntredR}(cnt));
                    else
                        tiling{cnt-1}{2,cnt2}.lowpass = gdct_cos((tiling{cnt-1}{2,cnt2}.polar_a-(cnt2-1+cntredA/red(1))*angles{1,cntredR}(cnt))/angles{1,cntredR}(cnt));
                    end
                    tiling{cnt-1}{2,cnt2}.hipass = sqrt(1-(tiling{cnt-1}{2,cnt2}.lowpass).^2);
                    tiling{cnt-1}{2,cnt2}.lowpass = (tiling{cnt-1}{2,cnt2}.lowpass).*(level{cnt-1}.hipass(pos));
                    tiling{cnt-1}{2,cnt2}.hipass = (tiling{cnt-1}{2,cnt2}.hipass).*(level{cnt-1}.hipass(pos));
                else
                    pos = find(level{cnt}.polar_a>=(cnt2-1+cntredA/red(1))*angles{1,cntredR}(cnt) | level{cnt}.polar_a < (cnt2+cntredA/red(1))*angles{1,cntredR}(cnt)-2*pi);
                    tiling{cnt-1}{1,cnt2} = struct('polar_a',0,'index_1',0,'index_2',0,'X_val',0,'X1_val',0,'X2_val',0,'lowpass',0,'hipass',0);
                    tiling{cnt-1}{1,cnt2}.polar_a = level{cnt}.polar_a(pos);
                    tiling{cnt-1}{1,cnt2}.index_1 = level{cnt}.index_1(pos);
                    tiling{cnt-1}{1,cnt2}.index_2 = level{cnt}.index_2(pos);
                    tiling{cnt-1}{1,cnt2}.X_val = level{cnt}.X_val(pos);
                    tiling{cnt-1}{1,cnt2}.X1_val = level{cnt}.X1_val(pos);
                    tiling{cnt-1}{1,cnt2}.X2_val = level{cnt}.X2_val(pos);
                    if ~is_cos
                        tiling{cnt-1}{1,cnt2}.lowpass = fdct_wrapping_window(mod(tiling{cnt-1}{1,cnt2}.polar_a-(cnt2-1+cntredA/red(1))*angles{1,cntredR}(cnt),2*pi)/angles{1,cntredR}(cnt));
                    else
                        tiling{cnt-1}{1,cnt2}.lowpass = gdct_cos(mod(tiling{cnt-1}{1,cnt2}.polar_a-(cnt2-1+cntredA/red(1))*angles{1,cntredR}(cnt),2*pi)/angles{1,cntredR}(cnt));
                    end
                    tiling{cnt-1}{1,cnt2}.hipass = sqrt(1-(tiling{cnt-1}{1,cnt2}.lowpass).^2);
                    tiling{cnt-1}{1,cnt2}.lowpass = (tiling{cnt-1}{1,cnt2}.lowpass).*(level{cnt}.lowpass(pos));
                    tiling{cnt-1}{1,cnt2}.hipass = (tiling{cnt-1}{1,cnt2}.hipass).*(level{cnt}.lowpass(pos));
                    pos = find(level{cnt-1}.polar_a>=(cnt2-1+cntredA/red(1))*angles{1,cntredR}(cnt) | level{cnt-1}.polar_a<(cnt2+cntredA/red(1))*angles{1,cntredR}(cnt)-2*pi);
                    tiling{cnt-1}{2,cnt2} = struct('polar_a',0,'index_1',0,'index_2',0,'X_val',0,'X1_val',0,'X2_val',0,'lowpass',0,'hipass',0);
                    tiling{cnt-1}{2,cnt2}.polar_a = level{cnt-1}.polar_a(pos);
                    tiling{cnt-1}{2,cnt2}.index_1 = level{cnt-1}.index_1(pos);
                    tiling{cnt-1}{2,cnt2}.index_2 = level{cnt-1}.index_2(pos);
                    tiling{cnt-1}{2,cnt2}.X_val = level{cnt-1}.X_val(pos);
                    tiling{cnt-1}{2,cnt2}.X1_val = level{cnt-1}.X1_val(pos);
                    tiling{cnt-1}{2,cnt2}.X2_val = level{cnt-1}.X2_val(pos);
                    if ~is_cos
                        tiling{cnt-1}{2,cnt2}.lowpass = fdct_wrapping_window(mod(tiling{cnt-1}{2,cnt2}.polar_a-(cnt2-1+cntredA/red(1))*angles{1,cntredR}(cnt),2*pi)/angles{1,cntredR}(cnt));
                    else
                        tiling{cnt-1}{2,cnt2}.lowpass = gdct_cos(mod(tiling{cnt-1}{2,cnt2}.polar_a-(cnt2-1+cntredA/red(1))*angles{1,cntredR}(cnt),2*pi)/angles{1,cntredR}(cnt));
                    end
                    tiling{cnt-1}{2,cnt2}.hipass = sqrt(1-(tiling{cnt-1}{2,cnt2}.lowpass).^2);
                    tiling{cnt-1}{2,cnt2}.lowpass = (tiling{cnt-1}{2,cnt2}.lowpass).*(level{cnt-1}.hipass(pos));
                    tiling{cnt-1}{2,cnt2}.hipass = (tiling{cnt-1}{2,cnt2}.hipass).*(level{cnt-1}.hipass(pos));
                end
            end
        end
        tiling{1} = cell(1,1);
        %maybe can speed up here!!!
        %generate C
        if R_low == 0
            temp = zeros(sz);
            temp((mod(level{1}.index_1,sz(1))+1)+((mod(level{1}.index_2,sz(2))))*sz(1)) = level{1}.X_val;
            temp((mod(level{2}.index_1,sz(1))+1)+((mod(level{2}.index_2,sz(2))))*sz(2)) = (level{2}.X_val).*(level{2}.lowpass);
            fac = 1;%sqrt(pi*wedge_end{1,cntredR}(2)^2);
            C = fftshift(ifft2(ifftshift(temp)))*sqrt(sz(1)*sz(2));
            temp = zeros(sz);
            temp((mod(level{1}.index_1,sz(1))+1)+((mod(level{1}.index_2,sz(2))))*sz(1)) = level{1}.X1_val;
            temp((mod(level{2}.index_1,sz(1))+1)+((mod(level{2}.index_2,sz(2))))*sz(2)) = (level{2}.X1_val).*(level{2}.lowpass);
            C1 = fftshift(ifft2(ifftshift(temp)))*sqrt(sz(1)*sz(2));
            temp = zeros(sz);
            temp((mod(level{1}.index_1,sz(1))+1)+((mod(level{1}.index_2,sz(2))))*sz(1)) = level{1}.X2_val;
            temp((mod(level{2}.index_1,sz(1))+1)+((mod(level{2}.index_2,sz(2))))*sz(2)) = (level{2}.X2_val).*(level{2}.lowpass);
            C2 = fftshift(ifft2(ifftshift(temp)))*sqrt(sz(1)*sz(2));
            % kk1 local wave vector in the first direction
            % kk2 local wave vector in the second direction
            gud = find(abs(C)>epsl);
            C1 = real( C1(gud)./C(gud) / (2*pi*i));
            C2 = real( C2(gud)./C(gud) / (2*pi*i));
            dist = sqrt(C1.^2+C2.^2);
            good = find(dist<=min(N/sqrt(2),R_high) & dist>=R_low);
            
            kk1 = repmat(EXT,sz);
            kk1(gud(good)) = C1(good);
            kk2 = kk1;
            kk2(gud(good)) = C2(good);
            
            for a=1:sz(1)
                for b=1:sz(2)
                    tc = C(a,b);
                    tk1 = kk1(a,b);
                    tk2 = kk2(a,b);
                    tc = tc(:);
                    gud = find(tk1<EXT);
                    tk1 = tk1(gud);
                    tk2 = tk2(gud);
                    for g=1:length(gud)
                        loc1 = round(tk1(g)*WB)+NB+1;
                        loc2 = round(tk2(g)*WB)+NB+1;
                        temp_energy = abs(tc(gud(g))).^2;
                        kb(loc1,loc2,a,b) = kb(loc1,loc2,a,b) + temp_energy;
                        avgdx(loc1,loc2,a,b) = avgdx(loc1,loc2,a,b) + tk1(g) * temp_energy;
                        avgdy(loc1,loc2,a,b) = avgdy(loc1,loc2,a,b) + tk2(g) * temp_energy;
                    end
                end
            end
        end
        for cnt =3:nbscales{1,cntredR}
            fac = 1;%sqrt(angles{1,cntredR}(cnt)*(wedge_end{1,cntredR}(cnt)^2-wedge_end{1,cntredR}(cnt-2)^2)/4);
            if angleRange(1)+2*pi-angleRange(2)<2*angles{1,cntredR}(cnt)
                cntRange = 1:nbangles{1,cntredR}(cnt);
            else
                cnt_st = mod(floor(angleRange(1)/angles{1,cntredR}(cnt)-cntredA/red(1))-1,nbangles{1,cntredR}(cnt))+1;
                cnt_ed = mod(ceil(angleRange(2)/angles{1,cntredR}(cnt)-cntredA/red(1))-1,nbangles{1,cntredR}(cnt))+1;
                if cnt_st<=cnt_ed
                    cntRange = cnt_st:cnt_ed;
                else
                    cntRange = [1:cnt_ed cnt_st:nbangles{1,cntredR}(cnt)];
                end
            end
            for cnt2 = cntRange
                %for cnt2 = 1:nbangles{1,cntredR}(cnt)
                temp = zeros(sz);
                x1 = mod(tiling{cnt-1}{1,cnt2}.index_1+pN1+1,sz(1))+1;
                x2 = mod(tiling{cnt-1}{1,cnt2}.index_2+pN2+1,sz(2))+1;
                
                %tmpf(tiling{cnt-1}{1,cnt2}.index_1+pN1+1+(tiling{cnt-1}{1,cnt2}.index_2+pN2+1-1)*N1) = 1;
                
                temp(x1+(x2-1)*sz(1)) = (tiling{cnt-1}{1,cnt2}.lowpass).*(tiling{cnt-1}{1,cnt2}.X_val);
                x1 = mod(tiling{cnt-1}{2,cnt2}.index_1+pN1+1,sz(1))+1;
                x2 = mod(tiling{cnt-1}{2,cnt2}.index_2+pN2+1,sz(2))+1;
                
                %tmpf(tiling{cnt-1}{2,cnt2}.index_1+pN1+1+(tiling{cnt-1}{2,cnt2}.index_2+pN2+1-1)*N1) = 1;
                
                temp(x1+(x2-1)*sz(1)) = (tiling{cnt-1}{2,cnt2}.lowpass).*(tiling{cnt-1}{2,cnt2}.X_val);
                if cnt2 < nbangles{1,cntredR}(cnt)
                    x1 = mod(tiling{cnt-1}{1,cnt2+1}.index_1+pN1+1,sz(1))+1;
                    x2 = mod(tiling{cnt-1}{1,cnt2+1}.index_2+pN2+1,sz(2))+1;
                    
                    %tmpf(tiling{cnt-1}{1,cnt2}.index_1+pN1+1+(tiling{cnt-1}{1,cnt2}.index_2+pN2+1-1)*N1) = 1;
                    
                    temp(x1+(x2-1)*sz(1)) = (tiling{cnt-1}{1,cnt2+1}.hipass).*(tiling{cnt-1}{1,cnt2+1}.X_val);
                    x1 = mod(tiling{cnt-1}{2,cnt2+1}.index_1+pN1+1,sz(1))+1;
                    x2 = mod(tiling{cnt-1}{2,cnt2+1}.index_2+pN2+1,sz(2))+1;
                    
                    %tmpf(tiling{cnt-1}{2,cnt2}.index_1+pN1+1+(tiling{cnt-1}{2,cnt2}.index_2+pN2+1-1)*N1) = 1;
                    
                    temp(x1+(x2-1)*sz(1)) = (tiling{cnt-1}{2,cnt2+1}.hipass).*(tiling{cnt-1}{2,cnt2+1}.X_val);
                else
                    x1 = mod(tiling{cnt-1}{1,1}.index_1+pN1+1,sz(1))+1;
                    x2 = mod(tiling{cnt-1}{1,1}.index_2+pN2+1,sz(2))+1;
                    
                    %tmpf(tiling{cnt-1}{1,cnt2}.index_1+pN1+1+(tiling{cnt-1}{1,cnt2}.index_2+pN2+1-1)*N1) = 1;
                    
                    temp(x1+(x2-1)*sz(1)) = (tiling{cnt-1}{1,1}.hipass).*(tiling{cnt-1}{1,1}.X_val);
                    x1 = mod(tiling{cnt-1}{2,1}.index_1+pN1+1,sz(1))+1;
                    x2 = mod(tiling{cnt-1}{2,1}.index_2+pN2+1,sz(2))+1;
                    
                    %tmpf(tiling{cnt-1}{2,cnt2}.index_1+pN1+1+(tiling{cnt-1}{2,cnt2}.index_2+pN2+1-1)*N1) = 1;
                    
                    temp(x1+(x2-1)*sz(1)) = (tiling{cnt-1}{2,1}.hipass).*(tiling{cnt-1}{2,1}.X_val);
                end
                C = fftshift(ifft2(ifftshift(temp)))*sqrt(sz(1)*sz(2));
                
                
                temp = zeros(sz);
                x1 = mod(tiling{cnt-1}{1,cnt2}.index_1+pN1+1,sz(1))+1;
                x2 = mod(tiling{cnt-1}{1,cnt2}.index_2+pN2+1,sz(2))+1;
                temp(x1+(x2-1)*sz(1)) = (tiling{cnt-1}{1,cnt2}.lowpass).*(tiling{cnt-1}{1,cnt2}.X1_val);
                x1 = mod(tiling{cnt-1}{2,cnt2}.index_1+pN1+1,sz(1))+1;
                x2 = mod(tiling{cnt-1}{2,cnt2}.index_2+pN2+1,sz(2))+1;
                temp(x1+(x2-1)*sz(1)) = (tiling{cnt-1}{2,cnt2}.lowpass).*(tiling{cnt-1}{2,cnt2}.X1_val);
                if cnt2 < nbangles{1,cntredR}(cnt)
                    x1 = mod(tiling{cnt-1}{1,cnt2+1}.index_1+pN1+1,sz(1))+1;
                    x2 = mod(tiling{cnt-1}{1,cnt2+1}.index_2+pN2+1,sz(2))+1;
                    temp(x1+(x2-1)*sz(1)) = (tiling{cnt-1}{1,cnt2+1}.hipass).*(tiling{cnt-1}{1,cnt2+1}.X1_val);
                    x1 = mod(tiling{cnt-1}{2,cnt2+1}.index_1+pN1+1,sz(1))+1;
                    x2 = mod(tiling{cnt-1}{2,cnt2+1}.index_2+pN2+1,sz(2))+1;
                    temp(x1+(x2-1)*sz(1)) = (tiling{cnt-1}{2,cnt2+1}.hipass).*(tiling{cnt-1}{2,cnt2+1}.X1_val);
                else
                    x1 = mod(tiling{cnt-1}{1,1}.index_1+pN1+1,sz(1))+1;
                    x2 = mod(tiling{cnt-1}{1,1}.index_2+pN2+1,sz(2))+1;
                    temp(x1+(x2-1)*sz(1)) = (tiling{cnt-1}{1,1}.hipass).*(tiling{cnt-1}{1,1}.X1_val);
                    x1 = mod(tiling{cnt-1}{2,1}.index_1+pN1+1,sz(1))+1;
                    x2 = mod(tiling{cnt-1}{2,1}.index_2+pN2+1,sz(2))+1;
                    temp(x1+(x2-1)*sz(1)) = (tiling{cnt-1}{2,1}.hipass).*(tiling{cnt-1}{2,1}.X1_val);
                end
                C1 = fftshift(ifft2(ifftshift(temp)))*sqrt(sz(1)*sz(2));
                temp = zeros(sz);
                x1 = mod(tiling{cnt-1}{1,cnt2}.index_1+pN1+1,sz(1))+1;
                x2 = mod(tiling{cnt-1}{1,cnt2}.index_2+pN2+1,sz(2))+1;
                temp(x1+(x2-1)*sz(1)) = (tiling{cnt-1}{1,cnt2}.lowpass).*(tiling{cnt-1}{1,cnt2}.X2_val);
                x1 = mod(tiling{cnt-1}{2,cnt2}.index_1+pN1+1,sz(1))+1;
                x2 = mod(tiling{cnt-1}{2,cnt2}.index_2+pN2+1,sz(2))+1;
                temp(x1+(x2-1)*sz(1)) = (tiling{cnt-1}{2,cnt2}.lowpass).*(tiling{cnt-1}{2,cnt2}.X2_val);
                if cnt2 < nbangles{1,cntredR}(cnt)
                    x1 = mod(tiling{cnt-1}{1,cnt2+1}.index_1+pN1+1,sz(1))+1;
                    x2 = mod(tiling{cnt-1}{1,cnt2+1}.index_2+pN2+1,sz(2))+1;
                    temp(x1+(x2-1)*sz(1)) = (tiling{cnt-1}{1,cnt2+1}.hipass).*(tiling{cnt-1}{1,cnt2+1}.X2_val);
                    x1 = mod(tiling{cnt-1}{2,cnt2+1}.index_1+pN1+1,sz(1))+1;
                    x2 = mod(tiling{cnt-1}{2,cnt2+1}.index_2+pN2+1,sz(2))+1;
                    temp(x1+(x2-1)*sz(1)) = (tiling{cnt-1}{2,cnt2+1}.hipass).*(tiling{cnt-1}{2,cnt2+1}.X2_val);
                else
                    x1 = mod(tiling{cnt-1}{1,1}.index_1+pN1+1,sz(1))+1;
                    x2 = mod(tiling{cnt-1}{1,1}.index_2+pN2+1,sz(2))+1;
                    temp(x1+(x2-1)*sz(1)) = (tiling{cnt-1}{1,1}.hipass).*(tiling{cnt-1}{1,1}.X2_val);
                    x1 = mod(tiling{cnt-1}{2,1}.index_1+pN1+1,sz(1))+1;
                    x2 = mod(tiling{cnt-1}{2,1}.index_2+pN2+1,sz(2))+1;
                    temp(x1+(x2-1)*sz(1)) = (tiling{cnt-1}{2,1}.hipass).*(tiling{cnt-1}{2,1}.X2_val);
                end
                C2 = fftshift(ifft2(ifftshift(temp)))*sqrt(sz(1)*sz(2));
                gud = find(abs(C)>epsl);
                C1 = real( C1(gud)./C(gud) / (2*pi*i));
                C2 = real( C2(gud)./C(gud) / (2*pi*i));
                dist = sqrt(C1.^2+C2.^2);
                good = find(dist<=min(N/sqrt(2),R_high) & dist>=R_low);
                
                kk1 = repmat(EXT,sz);
                kk1(gud(good)) = C1(good);
                kk2 = kk1;
                kk2(gud(good)) = C2(good);
                
                for a=1:sz(1)
                    for b=1:sz(2)
                        tc = C(a,b);
                        tk1 = kk1(a,b);
                        tk2 = kk2(a,b);
                        tc = tc(:);
                        gud = find(tk1<EXT);
                        tk1 = tk1(gud);
                        tk2 = tk2(gud);
                        for g=1:length(gud)
                            loc1 = round(tk1(g)*WB)+NB+1;
                            loc2 = round(tk2(g)*WB)+NB+1;
                            temp_energy = abs(tc(gud(g))).^2;
                            kb(loc1,loc2,a,b) = kb(loc1,loc2,a,b) + temp_energy;
                            avgdx(loc1,loc2,a,b) = avgdx(loc1,loc2,a,b) + tk1(g) * temp_energy;
                            avgdy(loc1,loc2,a,b) = avgdy(loc1,loc2,a,b) + tk2(g) * temp_energy;
                        end
                    end
                end
            end
        end
    end
end


%imagesc(tmpf);axis image;



