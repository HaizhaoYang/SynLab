function [realsize R_lowh R_highh R_lowv R_highv] = imgSize(N,M,name)
% This code stores the frequency range for SST
%
% By Haizhao Yang

switch name
    case 'synthetic'
        realsize = [1 1];
        R_lowh = 200; R_highh = 370;
        R_lowv = 200; R_highv = 370;
    case 'F205'
        realsize = [N*2.54/600 M*2.54/600];
        R_lowh = 180; R_highh = 270;
        R_lowv = 200; R_highv = 370;
    case 'F260'
        realsize = [N*2.54/600 M*2.54/600];
        R_lowh = 150; R_highh = 300;
        R_lowv = 200; R_highv = 350;
    case 'F458'
        realsize = [N*2.54/600 M*2.54/600];
        R_lowh = 250; R_highh = 400;
        R_lowv = 150; R_highv = 250;
    case 'F634'
        realsize = [N*2.54/600 M*2.54/600];
        R_lowh = 240; R_highh = 420;
        R_lowv = 160; R_highv = 240;
    case'F651'
        realsize = [N*2.54/600 M*2.54/600];
        R_lowh = 200; R_highh = 350;
        R_lowv = 150; R_highv = 250;
    case 'F659'
        realsize = [N*2.54/600 M*2.54/600];
        R_lowh = 250; R_highh = 340;
        R_lowv = 156; R_highv = 244;
    case 'L07'
        realsize = [N*2.54/600 M*2.54/600];
        R_lowh = 166; R_highh = 370;
        R_lowv = 166; R_highv = 370;
    case 'L11'
        realsize = [N*2.54/600 M*2.54/600];
        R_lowh = 160; R_highh = 320;
        R_lowv = 160; R_highv = 320;
    case 'L17'
        realsize = [N*2.54/600 M*2.54/600];
        R_lowh = 180; R_highh = 350;
        R_lowv = 180; R_highv = 350;
    case 'L30'
        realsize = [N*2.54/600 M*2.54/600];
        R_lowh = 200; R_highh = 350;
        R_lowv = 160; R_highv = 320;
    case 'D_Imitator'
        realsize = [N*82.1/8800,N*65.9/7080];
        R_lowh = 230; R_highh = 370;
        R_lowv = 230; R_highv = 370;
    case 'D_Rembrandt'
        realsize = [N*118.11/8800,N*96.52/6180];
        R_lowh = 170; R_highh = 340;
        R_lowv = 170; R_highv = 340;
    case 'D_Roussea'
        realsize = [N*99.7/8800,N*80/7080];
        R_lowh = 160; R_highh = 650;
        R_lowv = 160; R_highv = 650;
    case 'D_Giottoa'
        realsize = [N*42.822876/8800,N*35.08375/7020];
        R_lowh = 160; R_highh = 450;%weft
        R_lowv = 200; R_highv = 430;%warp
    case 'D_Giottob'
        realsize = [N*42.822876/8800,N*35.08375/7020];
        R_lowh = 150; R_highh = 340;%weft
        R_lowv = 240; R_highv = 430;%warp
    case 'D_Giottoc'
        realsize = [N*62.07125/8304,N*43.815/6096];
        R_lowh = 220*2; R_highh = 430*2;%warp
        R_lowv = 150*2; R_highv = 340*2;%weft
    case 'D_Giottod'
        realsize = [N*42.875/8800,N*35.08375/7020];
        R_lowh = 150; R_highh = 410;%weft
        R_lowv = 240; R_highv = 430;%warp
    case 'D_Giottoe'
        realsize = [N*42.743374/8800,N*35.08375/7020];
        R_lowh = 150; R_highh = 340;%weft
        R_lowv = 240; R_highv = 430;%warp
    case 'test'
        realsize = [N*1/2048,N/2048];
        R_lowh = 100*2; R_highh = 200*2;%weft
        R_lowv = 150*2; R_highv = 250*2;%warp
    case 'D_Giottonewa'
        realsize = [N*2.54/300 N*2.54/300];
        R_lowh = 350; R_highh = 700;%weft
        R_lowv = 370; R_highv = 700;%warp
    case 'D_Giottonewb'
        realsize = [N*2.54/300 N*2.54/300];
        R_lowh = 270; R_highh = 500;%weft
        R_lowv = 370; R_highv = 750;%warp
    case 'D_Giottonewc'
        realsize = [N*2.54/300 N*2.54/300];
        R_lowh = 370; R_highh = 700;%warp
        R_lowv = 300; R_highv = 600;%weft
    case 'D_Giottonewd'
        realsize = [N*2.54/300 N*2.54/300];
        R_lowh = 290; R_highh = 550;%weft
        R_lowv = 370; R_highv = 750;%warp
    case 'D_Giottonewe'
        realsize = [N*2.54/300 N*2.54/300];
        R_lowh = 290; R_highh = 550;%weft
        R_lowv = 370; R_highv = 770;%warp
    case 'Ryder'
        realsize = [N*2.54/2800,N*2.54/2800];
        R_lowh = 116; R_highh = 172;
        R_lowv = 110; R_highv = 156;
    case 'RyderPart'
        realsize = [N*2.54/2800,N*2.54/2800];
        R_lowh = 108; R_highh = 210;
        R_lowv = 102; R_highv = 192;
end
