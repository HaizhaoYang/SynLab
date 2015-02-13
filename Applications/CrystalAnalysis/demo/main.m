%% the code reads in the result from the synchrosqueezed wave packet transform applied to an image of a hexagonal crystal
%
% By Haizhao Yang

%% -------------------------------------------------------------------------
% Inputs:
% example       integer from 1 to 8
%
% Outputs:
% phi           original image
% G             optimized gradient field
% curlG         curl of G
% G0            initial gradient field provided by SST
% binaryIndex   fine defect regions with indices identifying connected
%               components

% mex LocBin.c;
% mex LocWeight.c;
% mex SS_polar.c;
close all;
is_Bin = 0;
is_ns = 1;
address = sprintf('results/Noise');
NM = 0.5;%-11.25 1.4; -2.35 0.5
%% -------------------------------------------------------------------------
% Generate data and choose parameters
%
% 1) parameters for frequency band identification
% 1.1) 'extendRatio', usually takes 0-0.1, will influence the radii
%      estimation on the boundarie
% 1.2) 'threRatio', a parameter to combine several possible close energy
%      peaks; threRatio=0.1 by default
% 1.3) 'epsFreqBand', denoising threshold when identifying the frequency
%      band; epsFreqBand=1e-2 by default
%
% 2) parameters for the synchrosqueezed transform
% 2.1) 'numWave', the number of local wave vectors of interest, determined
%      by the crystal, e.g., numWave=3 for hexagonal crystals
% 2.2) 'numGridPerSector', the number of grid points in one sector;
%      in the case of hexagonal crystals, each sector covers an angle of
%      pi/3; numGridPerSector = 30 by default
% 2.3) 'red', a redundancy parameter of the wave packet transform; red=1 by
%      default
% 2.4) 'rad' control the size of the support of the mother wave packet in
%      the K-plane; rad in [0.5,2]
% 2.5) 's_sc' and 't_sc' are scaling parameters for wave packet transforms
% 2.6) 'epsSST', denoising threshold; epsSST=1e-2 by default
%
% 3) parameters for defect region identification
% 3.1) 'energyThre', a thresholding parameter to identify valleys from
%      energy landmass provided by SST; energyThre in [0,3], smaller vthre
%      leads to more narrow boundaries;
% 3.2) 'coarseThre', a thresholding parameter to generate a coarse binary
%      image identifying grains and grain boundaries;
%      coarseThre in [0,1]; smaller value leads to wider boundaries
% 3.3) 'fineThre', a thresholding parameter to generate a fine binary
%      image identifying grains and grain boundaries;
%      fineThre in [0,1]; larger value leads to more narrow boundaries
%
% 4) parameters for the optimization
% 4.1) 'numIter', number of interations;
% 4.2) 'disVolThre', a thresholding parameter to show the distortion
%      volumne;
example = 1;
switch example
    case 1 % use high subsampleRate or small red
        close all;
        load 'data/GB1.mat';
        if is_ns
            if NM>1
                parameters = struct('extendRatio',0.0,'threRatio',0.1, ...
                    'epsFreqBand',1e-2,'numWave',3, ...
                    'numGridPerSector',30,'red',8,'rad',2, ...
                    's_sc',0.6,'t_sc',0.6,'epsSST',0.15,'is_Bin',is_Bin);
                subsampleRate = 4;
            else if NM>0.5
                    parameters = struct('extendRatio',0.0,'threRatio',0.1, ...
                        'epsFreqBand',1e-2,'numWave',3, ...
                        'numGridPerSector',30,'red',10,'rad',2, ...
                        's_sc',0.625,'t_sc',0.625,'epsSST',0.15,'is_Bin',is_Bin);
                    subsampleRate = 4;% NM <= 1
                else
                    parameters = struct('extendRatio',0.0,'threRatio',0.1, ...
                        'epsFreqBand',1e-2,'numWave',3, ...
                        'numGridPerSector',30,'red',10,'rad',2, ...
                        's_sc',0.625,'t_sc',0.625,'epsSST',0.15,'is_Bin',is_Bin);
                    subsampleRate = 4;% NM <= 1
                end
            end
        else
            parameters = struct('extendRatio',0.0,'threRatio',0.1, ...
                'epsFreqBand',1e-2,'numWave',3, ...
                'numGridPerSector',30,'red',4,'rad',2, ...
                's_sc',0.8,'t_sc',0.8,'epsSST',0.1,'is_Bin',is_Bin);
            subsampleRate = 2;
        end
        
        if ~is_ns
            energyThre = 2.3;%clean
        else
            energyThre = 2;%noisy
        end
        disVolThre = 0.1;
        numIter = 40;
end

% normalize the given image so as to use uniform thresholding parameters
fff = phi;
fff = fff/max(abs(fff(:)));
fff = fff - mean(fff(:));
if is_ns
    [Nx,Ny] = size(fff);
    ns = real(NM*(randn(Nx,Ny)+i*randn(Nx,Ny)));
    snr = 10*log10((var(fff(:)))/var(ns(:)))
    fff = fff + ns;
end

%% -------------------------------------------------------------------------
% The synchrosqueezed wave packet transform provids an initial guess
tic;
plotFIG = false;
[radii angles TTEng_1st TTEng_2nd masses] = initialGuess(fff,subsampleRate,plotFIG,parameters); % assumes y-coordinate pointing into direction of growing first index
BD = findGB(TTEng_1st,TTEng_2nd);

% bring wave vectors into order according to angle
[angles,perm] = sort(angles,3);
for k = 1:size(perm,1)
    for l = 1:size(perm,2)
        radii(k,l,:) = radii(k,l,perm(k,l,:));
    end
end

v = smoothImage( sum(masses,3), 10, 2 );
pos = find(v>energyThre);
N = size(BD,1);
binaryIndex = zeros(N,N);
binaryIndex(pos) = 1;

% find initial noisy deformation gradient, assuming y-coordinate pointing into direction of decreasing 1st index, and using two different point group interpretations
G = initGFromSSPeaks( angles, radii, pi/3*[0;1;2], median(radii(:))*ones(3,1) );

%% -------------------------------------------------------------------------
% show results
for ii = 1:2
    pic = figure;
    imagesc([0,1],[0,1],G(:,:,1,ii)); axis image;axis xy;
    %freezeColors;
%     colorbar;%cbfreeze(colorbar);
%     name = sprintf('%sG1%d_Ex%d.fig',address,ii,example);
%     saveas(pic,name);
end
for ii = 1:2
    pic = figure;
    imagesc([0,1],[0,1],G(:,:,2,ii)); axis image;axis xy;
    %freezeColors;
    colorbar;%cbfreeze(colorbar);
%     name = sprintf('%sG2%d_Ex%d.fig',address,ii,example);
%     saveas(pic,name);
end


% volume distortion
pic = figure;
vol = G(:,:,1).*G(:,:,4)-G(:,:,2).*G(:,:,3);
imagesc([0,1],[0,1],vol-1); axis image;  h=colorbar; caxis([-disVolThre disVolThre]);axis xy;
colorbar;
title('Distortion Volume');
% name = sprintf('%sVol_Ex%d.fig',address,example);
% saveas(pic,name);

mv = -92;
mv2 = 500;
% input image
pic = figure; hold on;
%1024-(93:156)-mv
imagesc([0,1],[0,1],fff);
imagesc([0,0.4],[0,0.4],real(fff(1024-163:1024-100,1:64)));
colorbar;
rectangle('Position',[(93+mv)/1024,861/1024,(256-193+1)/1024,(256-193+1)/1024]);
rectangle('Position',[0,0,0.4,0.4]);
hold off;
axis image;axis xy;
% name = sprintf('%sImg_Ex%d.fig',address,example);
% saveas(pic,name);

% grain boundary
[m n] = size(fff);
[mm nn] = size(BD);
rt = mm/m;
pic = figure; hold on;
imagesc([0,1],[0,1],BD);axis image;colorbar;axis xy; colormap (1-gray);
imagesc([0,0.4],[0,0.4],real(BD(ceil((1024-163:1024-100)*rt),ceil((1:64)*rt))));
colorbar;
r = rectangle('Position',[(mv+93)/1024,round(861)/1024,round((256-193+1))/1024,round((256-193+1))/1024]);
set(r,'edgecolor','k')
r = rectangle('Position',[0,0,0.4,0.4]);
set(r,'edgecolor','k')
hold off;
axis image;axis xy;
% name = sprintf('%sBD_Ex%d.fig',address,example);
% saveas(pic,name);
