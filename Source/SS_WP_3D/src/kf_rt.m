function r = kf_rt(w,n)
% kf_rt.m - right bump
%
% Written by Lexing Ying and Laurent Demanet, 2006
  
  r = zeros(size(w));
  an = pi/2*(n+1/2);
  en = (-1)^n;
  en1 = (-1)^(n+1);
  
  r = exp(-i*w/2) .* ( exp(i*an)*g(en*(w-pi*(n+1/2))) );

  
