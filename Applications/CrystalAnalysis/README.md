# CrystalTool

1. OVERVIEW

CrystalTool is a MATLAB toolbox to analyze atomic crystal images. For a given atomic crystal image, it contains several tools to identify grain boundary, crystal orientation, elastic deformation. A few examples of synthetic and real atomic crystal images are provided to illustrate how to use these tools.

2. INTRODUCTION

SSTmethod: A collection of MATLAB and MEX routines which applies 2D synchrosqueezed transforms proposed in [1][2] to analyze crystal images. Detailed description of this method is in [4]. An instruction for parameter tuning in the synchrosqueezed transform is in [3]. 

VarSSTmethod: A collection of MATLAB and MEX routines which applies a variational model to optimize the crystal analysis results provided by 2D synchrosqueezed transforms. Detailed description of this method is in [5]. 


[1] H. Yang and L. Ying. Synchrosqueezed wave packet transform for 2d mode decompo- sition. SIAM Journal on Imaging Sciences, 6(4):1979–2009, 2013.

[2] H. Yang and L. Ying. Synchrosqueezed curvelet transform for two-dimensional mode decomposition. SIAM Journal on Mathematical Analysis, 46(3):2052–2083, 2014.

[3] H. Yang and L. Ying. Robustness Analysis of Synchrosqueezed Transforms. Preprint.

[4] H. Yang, J. Lu, and L. Ying. Crystal image analysis using 2d synchrosqueezed trans- forms. arXiv:1402.1262 [math.NA], Submitted on 6 Feb 2014.

[5] J. Lu, B. Wirth and H. Yang. Compbining 2d synchrosqueezed wave packet transforms with optimization for crystal image analysis. arXiv:1501.06254 [cond-mat.mtrl-sci], Submitted on 26 Jan 2015.

3. INSTALLING SYNLAB

run the file SetPath.m first. It will automatically add all the MATLAB codes to your MATLAB path and compile all MEX files. After this, you can run all demo codes to see how to use this tool box.

4. COPY RIGHT

CrystalTool is copyright reserved. For further information, please contact 
Contact information
Jianfeng Lu at jianfeng@math.duke.edu
Benedikt Wirth at Benedikt.Wirth@uni-muenster.de
Haizhao Yang at haizhao@stanford.edu
Lexing Ying at lexing@math.stanford.edu