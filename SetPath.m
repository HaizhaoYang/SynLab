% Copyright (c) Haizhao Yang, Stanford University, 2014

global SetPath
global CSPT
global MSPT
		
type = computer;

if strcmp(type,'MAC2'),
  CSPT = ':';
  SetPath = [pwd, CSPT];
  MSPT = ';';
elseif isunix,
  % Mac OS X returns isunix=1
  CSPT = '/';
  SetPath = [pwd, CSPT];
  MSPT = ':';
elseif strcmp(type(1:2),'PC');
  CSPT = '\';	  
  SetPath = [pwd, CSPT];  
  MSPT = ';';
end

disp('Begin to set MATLAB path...')

back = CSPT;
tempPath = path;
front = [MSPT SetPath];
tempPath = [tempPath front];

% Source folder
tempPath = [tempPath front 'results' back];

tempPath = [tempPath front 'Source' back];

tempPath = [tempPath front 'Source' back 'SS_CT_2D' back];
tempPath = [tempPath front 'Source' back 'SS_CT_2D' back 'demo' back];
tempPath = [tempPath front 'Source' back 'SS_CT_2D' back 'src' back];

tempPath = [tempPath front 'Source' back 'SS_WP_3D' back];
tempPath = [tempPath front 'Source' back 'SS_WP_3D' back 'demo' back];
tempPath = [tempPath front 'Source' back 'SS_WP_3D' back 'src' back];

tempPath = [tempPath front 'Source' back 'SS_WP_1D' back];
tempPath = [tempPath front 'Source' back 'SS_WP_1D' back 'demo' back];
tempPath = [tempPath front 'Source' back 'SS_WP_1D' back 'src' back];

% Applications folder
tempPath = [tempPath front 'Applications' back];

tempPath = [tempPath front 'Applications' back 'CanvasAnalysis' back];
tempPath = [tempPath front 'Applications' back 'CanvasAnalysis' back 'demo' back];
tempPath = [tempPath front 'Applications' back 'CanvasAnalysis' back 'src' back];

tempPath = [tempPath front 'Applications' back 'SynCrystal' back];
tempPath = [tempPath front 'Applications' back 'SynCrystal' back 'data' back];
tempPath = [tempPath front 'Applications' back 'SynCrystal' back 'SSTmethod' back];
tempPath = [tempPath front 'Applications' back 'SynCrystal' back 'SSTmethod' back 'demo' back];
tempPath = [tempPath front 'Applications' back 'SynCrystal' back 'SSTmethod' back 'src' back];
tempPath = [tempPath front 'Applications' back 'SynCrystal' back 'VarSSTmethod' back];
tempPath = [tempPath front 'Applications' back 'SynCrystal' back 'VarSSTmethod' back 'demo' back];
tempPath = [tempPath front 'Applications' back 'SynCrystal' back 'VarSSTmethod' back 'src' back];
tempPath = [tempPath front 'Applications' back 'SynCrystal' back 'VarSSTmethod' back 'src' back 'srcOptPart' back];
tempPath = [tempPath front 'Applications' back 'SynCrystal' back 'VarSSTmethod' back 'src' back 'srcOther' back];
tempPath = [tempPath front 'Applications' back 'SynCrystal' back 'VarSSTmethod' back 'src' back 'srcSST' back];

tempPath = [tempPath front 'Applications' back 'DemoRobustness' back];
tempPath = [tempPath front 'Applications' back 'DemoRobustness' back 'data' back];
tempPath = [tempPath front 'Applications' back 'DemoRobustness' back 'demo' back];

tempPath = [tempPath front 'Applications' back 'BenchMarkSignals' back];
tempPath = [tempPath front 'Applications' back 'BenchMarkSignals' back 'demo' back];
tempPath = [tempPath front 'Applications' back 'BenchMarkSignals' back 'data' back];

tempPath = [tempPath front 'Applications' back 'SeismicWaveField' back];
tempPath = [tempPath front 'Applications' back 'SeismicWaveField' back 'data' back];
tempPath = [tempPath front 'Applications' back 'SeismicWaveField' back 'src' back];
tempPath = [tempPath front 'Applications' back 'SeismicWaveField' back 'demo' back];

tempPath = [tempPath front 'Applications' back 'GeneralModeDecom1D' back];
tempPath = [tempPath front 'Applications' back 'GeneralModeDecom1D' back 'data' back];
tempPath = [tempPath front 'Applications' back 'GeneralModeDecom1D' back 'src' back];
tempPath = [tempPath front 'Applications' back 'GeneralModeDecom1D' back 'demo' back];

tempPath = [tempPath front 'Applications' back 'CanvasRemoval' back];
tempPath = [tempPath front 'Applications' back 'CanvasRemoval' back 'data' back];
tempPath = [tempPath front 'Applications' back 'CanvasRemoval' back 'src' back];

path(tempPath);

disp('Begin to compile MEX files...');
rootDir = pwd;
cd(['Source' CSPT 'SS_CT_2D' CSPT 'src' CSPT]);
mex SS_polar.c;
mex SS_polar_old.c;
mex SS_polar_v2.c;
mex SS_polar_v1.c;
cd(rootDir);
cd(['Applications' CSPT 'SynCrystal' CSPT 'SSTmethod' CSPT 'src' CSPT]);
mex LocWeight.c;
mex LocSmooth.c;
cd(rootDir);
cd(['Applications' CSPT 'SynCrystal' CSPT 'VarSSTmethod' CSPT 'src' CSPT 'srcSST' CSPT]);
mex LocWeight.c;
mex LocWeight_v2.c;
cd(rootDir);
cd(['Applications' CSPT 'CanvasAnalysis' CSPT 'src' CSPT]);
mex LocWavVec_v2.c;
cd(rootDir);

disp('Path set!');

clear tempPath front back
clear SetPath MATLABVERSION CSPT
clear type MSPT

