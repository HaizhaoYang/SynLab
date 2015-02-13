function r = g(w)
% g.m - 'g' function Villemoes's construction
% 
% Written by Lexing Ying and Laurent Demanet, 2006
  r = zeros(size(w));
  gd = w<5*pi/6 & w>-7*pi/6;
  r(gd) = abs(sf(w(gd)-3*pi/2));
  
%----------------------------------------------------------------------
function r = sf(w)
  r = zeros(size(w));
  
  aw = abs(w);
  r(aw<=2*pi/3) = 0;
  
  gd = aw>=2*pi/3 & aw<=4*pi/3;
  r(gd) = 1/sqrt(2) * hf(w(gd)/2+pi);
  
  gd = aw>=4*pi/3 & aw<=8*pi/3;
  r(gd) = 1/sqrt(2) * hf(w(gd)/4);
  
  r(aw>8*pi/2) = 0;
  
function r = hf(w)
  w = mod(w+pi,2*pi) - pi;
  r = zeros(size(w));
  w = abs(w);
  
  r = sqrt(2) * cos(pi/2 * beta(3*w/pi-1));
  r(w<=pi/3) = sqrt(2);
  r(w>=2*pi/3) = 0;
  
function r = beta(x)
  r = (1-cos(pi*x))/2;
  %t = (1-cos(pi*x))/2;  r = (1-cos(pi*t))/2;
  %r = x.^4 .*(35-84*x + 70*x.^2 -20*x.^3);
  
    
  
