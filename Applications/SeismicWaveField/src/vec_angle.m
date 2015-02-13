function angle = vec_angle(a1,a2,b1,b2)
%Return the angle between two vectors (a1,a2) and (b1,b2)
dot_rsl = a1.*b1+a2.*b2;
norm_rsl = sqrt(a1.^2+a2.^2).*norm([b1 b2]);
angle = acos(dot_rsl./norm_rsl);

