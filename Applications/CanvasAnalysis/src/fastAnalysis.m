function [HorThre,VerThre,HorAgl,VerAgl,Xwhole,orgSize] = fastAnalysis(name,t_sc,s_sc,epsl,rad,red,NL,numGrid,numdiv,numdiv2,div)
%==========================================================================
% This code extracts canvas texture from a real painting and uses 2D
% synchrosqueezed transforms to count threads in the canvas texture.
%
% Input:
%   name    name of the data document
%   t_sc    scaling parameter for radius in SST
%               [default set to 1-1/8]
%   s_sc    scaling parameter for angle in SST
%               [default set to 1/2+1/8]
%   epsl    a threshold in SST, [default set to 1e-4] 
%   rad     the smallest size of supports of wave packets
%               [default set to 1]
%   red     the redundancy parameter, red(1) for angle, red(2) for radii
%   numGrid, number of grid points in the angle domain of the SST
%   NL      controls the range of interesting thread angle, e.g., the range
%           of the thread angle is [-180*NL/numGrid, 180*NL/numGrid]
%           degrees
%   numdiv  a parameter for the size of a piece of canvas, if your RAM is
%           about 32GB, you could choose 2, otherwise, numdiv should be larger
%   numdiv2 a parameter to down sample the canvas to save memory, if your RAM is
%           about 32GB, you could choose 4, otherwise, numdiv2 should be larger
%   div     a parameter to tune the resolution of the SST in the spatial
%           domain, if your RAM is about 32GB, you could choose 4, otherwise, 
%           div should be larger
%
% Output:
%   HorThre  horizontal thread count
%   VerThre  vertical thread count
%   HorAgl   horizontal thread angle
%   VerAgl   vertical thread angle
%   Xwhole   the canvas image
%   orgSize  the size of the canvas image
%
% Remark: Usually, you could tune two parameters, [t_sc s_sc] according 
%         to the paper "Robustness Analysis of Synchrosqueezed Transforms",
%         H. Yang and L. Ying. For clean data, [t_sc s_sc] = [1 0.8]. For 
%         noisy data, [t_sc s_sc] = [0.8 0.625]. The analysis is not
%         sensitive to other parameters.
%
%
% By Haizhao Yang
%==========================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set up parameters
if nargin<2, t_sc = 1; end;
if nargin<3, s_sc = 0.8; end;
if nargin<4, epsl = 1e-4; end;
if nargin<5, rad = 1; end;
if nargin<6, red = [8 5]; end;
if nargin<7, NL = 15; end;
if nargin<8, numGrid = 90; end;
if nargin<9, numdiv = 2; end;
if nargin<10, numdiv2 = 4; end;
if nargin<11, div = 4; end;
is_real = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set up the image
fileName = sprintf('%s.mat',name);
load(fileName);

N = 4096/numdiv;

old = 0;
[realsize R_lowh R_highh R_lowv R_highv] = imgSize(N,N,name);

pic=figure;imagesc([0,1],[0,1],real(X));axis image;colorbar;title('Original image');xlabel('x_1');ylabel('x_2');colormap gray;
axis([0 1 0 1]);

R_lowh = R_lowh/numdiv; R_highh = R_highh/numdiv; R_lowv = R_lowv/numdiv; R_highv = R_highv/numdiv;

Xwhole = X;
[m,n] = size(Xwhole);
orgSize = [m,n];
X = zeros(ceil(m/(0.5*N))*0.5*N,ceil(n/(0.5*N))*0.5*N);
[mm nn] = size(X);

range1 = (1:m)+floor((mm-m)/2);
range2 = (1:n)+floor((nn-n)/2);
X(range1,range2) = Xwhole;
[rm rn] = size(X);
Xwhole = X;
[m,n] = size(Xwhole);
numm = m/0.5/N-1;
numn = n/0.5/N-1;
HorThreCell = cell(numm,numn);
VerThreCell = cell(numm,numn);
HorAglCell = cell(numm,numn);
VerAglCell = cell(numm,numn);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2D SST for thread counting
for cntm = 1:numm
    for cntn = 1:numn
        [cntm cntn]
        %set swatch
        X = Xwhole((1:N)+(cntm-1)*N/2,(1:N)+(cntn-1)*N/2);
        %down sampling the original image to reduce possible computational cost
        Y=fftshift(fft2(ifftshift(X)));
        Y = Y(N*(numdiv2-1)/2/numdiv2+1:end-N*(numdiv2-1)/2/numdiv2,N*(numdiv2-1)/2/numdiv2+1:end-N*(numdiv2-1)/2/numdiv2);
        X = ifftshift(ifft2(fftshift(Y)));
        X = X/max(max(abs(X)));
        
        %Compute horizontal thread count
        %set up size of samples in the space
        SPg = ceil(size(X)/div/8)*8;              %May be out of memory if use size(X)
        R_low = R_lowh; R_high = R_highh;
        NB = [(R_high - R_low),numGrid];
        num_direction = 1;
        
        angleRange = [-pi*NL/numGrid pi*NL/numGrid];
        [ss_energy avgdx avgdy] = SS_ct2_polar_v2(num_direction,X,SPg,NB,rad,is_real,[R_low,R_high],angleRange,epsl,red,t_sc, s_sc);
        SPg(1) = size(ss_energy,3);
        SPg(2) = size(ss_energy,4);
        SPgh = SPg;
        
        %count threads
        Thread_h = zeros(SPg);
        Angle_h = zeros(SPg);
        L = size(ss_energy,1);
        [Angle_h,Thread_h] = LocWavVec_v2(ss_energy,avgdx,avgdy,round(NB(2)/2+1-NL),round(NB(2)/2+NL));
        Angle_h = Angle_h-pi/2;
        xx = Thread_h.*cos(Angle_h)/realsize(1); yy = Thread_h.*sin(Angle_h)/realsize(2);
        Thread_h = sqrt(xx.^2+yy.^2);
        Angle_h = atan(yy./xx)*180/pi;

        clear ss_energy svgdx svgdy;
        
        %----------------------------------------------------------------------
        %Compute vertical thread count
        
        %set up parameters
        SPg = ceil(size(X)/div/8)*8;
        R_low = R_lowv; R_high = R_highv;
        NB = [(R_high - R_low),numGrid];
        angleRange = [pi/2-pi*NL/numGrid pi/2+pi*NL/numGrid];
        [ss_energy avgdx avgdy] = SS_ct2_polar_v2(num_direction,X,SPg,NB,rad,is_real,[R_low,R_high],angleRange,epsl,red,t_sc,s_sc);
        SPg(1) = size(ss_energy,3);
        SPg(2) = size(ss_energy,4);
        SPgv = SPg;
        
        %count threads
        Thread_v = zeros(SPg);
        Angle_v = zeros(SPg);
        L = size(ss_energy,1);
        [Angle_v Thread_v] = LocWavVec_v2(ss_energy,avgdx,avgdy,NB(2)-NL+1,NL);
        Angle_v = mod(Angle_v+pi/2,pi)-pi/2;
        xx = Thread_v.*cos(Angle_v)/realsize(2); yy = Thread_v.*sin(Angle_v)/realsize(1);
        Thread_v = sqrt(xx.^2+yy.^2);
        Angle_v = atan(yy./xx)*180/pi;

        clear ss_energy avgdx avgdy;
        
        HorThreCell{cntm,cntn} = Thread_h;
        VerThreCell{cntm,cntn} = Thread_v;
        HorAglCell{cntm,cntn} = Angle_h;
        VerAglCell{cntm,cntn} = Angle_v;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% generate results
HorThre = zeros(SPgh/2.*[numm+1 numn+1]);
VerThre = zeros(SPgv/2.*[numm+1 numn+1]);
HorAgl = zeros(SPgh/2.*[numm+1 numn+1]);
VerAgl = zeros(SPgv/2.*[numm+1 numn+1]);
for cntm = 1:numm
    for cntn = 1:numn
        if cntm == 1 & cntn == 1
            temp = HorThreCell{cntm,cntn};
            HorThre(1:SPgh(1),1:SPgh(2)) = temp(1:SPgh(1),1:SPgh(2));
            temp = VerThreCell{cntm,cntn};
            VerThre(1:SPgv(1),1:SPgv(2)) = temp(1:SPgv(1),1:SPgv(2));
            temp = HorAglCell{cntm,cntn};
            HorAgl(1:SPgh(1),1:SPgh(2)) = temp(1:SPgh(1),1:SPgh(2));
            temp = VerAglCell{cntm,cntn};
            VerAgl(1:SPgv(1),1:SPgv(2)) = temp(1:SPgv(1),1:SPgv(2));
        else if cntm ==1 & cntn ~= 1
                temp = HorThreCell{cntm,cntn};
                HorThre(1:SPgh(1),((SPgh(2)*0.25+1):SPgh(2))+0.5*SPgh(2)*(cntn-1)) = temp(1:SPgh(1),SPgh(2)*0.25+1:end);
                temp = VerThreCell{cntm,cntn};
                VerThre(1:SPgv(1),((SPgv(2)*0.25+1):SPgv(2))+0.5*SPgv(2)*(cntn-1)) = temp(1:SPgv(1),SPgv(2)*0.25+1:end);
                temp = HorAglCell{cntm,cntn};
                HorAgl(1:SPgh(1),((SPgh(2)*0.25+1):SPgh(2))+0.5*SPgh(2)*(cntn-1)) = temp(1:SPgh(1),SPgh(2)*0.25+1:end);
                temp = VerAglCell{cntm,cntn};
                VerAgl(1:SPgv(1),((SPgv(2)*0.25+1):SPgv(2))+0.5*SPgv(2)*(cntn-1)) = temp(1:SPgv(1),SPgv(2)*0.25+1:end);
            else if cntm ~= 1 & cntn == 1
                    temp = HorThreCell{cntm,cntn};
                    HorThre(((SPgh(1)*0.25+1):SPgh(1))+0.5*SPgh(1)*(cntm-1),1:SPgh(2)) = temp(SPgh(1)*0.25+1:end,1:SPgh(2));
                    temp = VerThreCell{cntm,cntn};
                    VerThre(((SPgv(1)*0.25+1):SPgv(1))+0.5*SPgv(1)*(cntm-1),1:SPgv(2)) = temp(SPgv(1)*0.25+1:end,1:SPgv(2));
                    temp = HorAglCell{cntm,cntn};
                    HorAgl(((SPgh(1)*0.25+1):SPgh(1))+0.5*SPgh(1)*(cntm-1),1:SPgh(2)) = temp(SPgh(1)*0.25+1:end,1:SPgh(2));
                    temp = VerAglCell{cntm,cntn};
                    VerAgl(((SPgv(1)*0.25+1):SPgv(1))+0.5*SPgv(1)*(cntm-1),1:SPgv(2)) = temp(SPgv(1)*0.25+1:end,1:SPgv(2));
                else
                    temp = HorThreCell{cntm,cntn};
                    HorThre(((SPgh(1)*0.25+1):SPgh(1))+0.5*SPgh(1)*(cntm-1),((SPgh(2)*0.25+1):SPgh(2))+0.5*SPgh(2)*(cntn-1)) = temp(SPgh(1)*0.25+1:end,SPgh(2)*0.25+1:end);
                    temp = VerThreCell{cntm,cntn};
                    VerThre(((SPgv(1)*0.25+1):SPgv(1))+0.5*SPgv(1)*(cntm-1),((SPgv(2)*0.25+1):SPgv(2))+0.5*SPgv(2)*(cntn-1)) = temp(SPgv(1)*0.25+1:end,SPgv(2)*0.25+1:end);
                    temp = HorAglCell{cntm,cntn};
                    HorAgl(((SPgh(1)*0.25+1):SPgh(1))+0.5*SPgh(1)*(cntm-1),((SPgh(2)*0.25+1):SPgh(2))+0.5*SPgh(2)*(cntn-1)) = temp(SPgh(1)*0.25+1:end,SPgh(2)*0.25+1:end);
                    temp = VerAglCell{cntm,cntn};
                    VerAgl(((SPgv(1)*0.25+1):SPgv(1))+0.5*SPgv(1)*(cntm-1),((SPgv(2)*0.25+1):SPgv(2))+0.5*SPgv(2)*(cntn-1)) = temp(SPgv(1)*0.25+1:end,SPgv(2)*0.25+1:end);
                end
            end
        end
    end
end
[m n] = size(HorThre);
HorThre = HorThre(  ceil(range1(1)*m/rm):floor(range1(end)*m/rm) , ceil(range2(1)*n/rn):floor(range2(end)*n/rn) );
HorAgl = HorAgl(  ceil(range1(1)*m/rm):floor(range1(end)*m/rm) , ceil(range2(1)*n/rn):floor(range2(end)*n/rn) );
[m,n] = size(VerThre);
VerThre = VerThre(  ceil(range1(1)*m/rm):floor(range1(end)*m/rm) , ceil(range2(1)*n/rn):floor(range2(end)*n/rn) );
VerAgl = VerAgl(  ceil(range1(1)*m/rm):floor(range1(end)*m/rm) , ceil(range2(1)*n/rn):floor(range2(end)*n/rn) );

pic=figure;
imagesc([0,1],[0,1],HorThre);axis image;colorbar;caxis([R_lowh R_highh]/realsize(1));head = sprintf('Horizontal Thread Count, mean=%f',mean(HorThre(:)));title(head);axis tight;
xlabel('x_1');ylabel('x_2');axis([0 1 0 1]); 

pic = figure;imagesc([0,1],[0,1],HorAgl);axis image;colorbar;head = sprintf('Horizontal Thread Angle, mean=%f',mean(HorAgl(:)));title(head);axis tight;
xlabel('x_1');ylabel('x_2');axis([0 1 0 1]); 

pic=figure;
imagesc([0,1],[0,1],VerThre);axis image;colorbar;caxis([R_lowv R_highv]/realsize(2));head = sprintf('Vertical Thread Count, mean=%f',mean(VerThre(:)));title(head);axis tight;
xlabel('x_1');ylabel('x_2');axis([0 1 0 1]); 

pic = figure;imagesc([0,1],[0,1],VerAgl);axis image;colorbar;head = sprintf('Vertical Thread Angle, mean=%f',mean(VerAgl(:)));title(head);axis tight;
xlabel('x_1');ylabel('x_2');axis([0 1 0 1]); 

HorThreDev = HorThre - mean(HorThre(:));
HorAglDev = HorAgl - mean(HorAgl(:));
VerThreDev = VerThre - mean(VerThre(:));
VerAglDev = VerAgl - mean(VerAgl(:));
pic=figure;
imagesc([0,1],[0,1],HorThreDev);axis image;colorbar;head = sprintf('Horizontal Thread Count Deviation, mean=%f',mean(HorThre(:)));title(head);axis tight;
xlabel('x_1');ylabel('x_2');axis([0 1 0 1]); 

pic = figure;imagesc([0,1],[0,1],HorAglDev);axis image;colorbar;head = sprintf('Horizontal Thread Angle Deviation, mean=%f',mean(HorAgl(:)));title(head);axis tight;
xlabel('x_1');ylabel('x_2');axis([0 1 0 1]); 

pic=figure;
imagesc([0,1],[0,1],VerThreDev);axis image;colorbar;head = sprintf('Vertical Thread Count Deviation, mean=%f',mean(VerThre(:)));title(head);axis tight;
xlabel('x_1');ylabel('x_2');axis([0 1 0 1]); 

pic = figure;imagesc([0,1],[0,1],VerAglDev);axis image;colorbar;head = sprintf('Vertical Thread Angle Deviation, mean=%f',mean(VerAgl(:)));title(head);axis tight;
xlabel('x_1');ylabel('x_2');axis([0 1 0 1]); 

save(sprintf('results/%sresults.mat',name));
