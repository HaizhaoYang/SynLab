function y = gen_shape2(x,num)
% Generate wave shapes

x = x - floor(x);
y = zeros(size(x));
if num == 1
    y = (cos(0.8*cos(x))-tan(x).*sin(0.8*cos(x))).*cos(x);
else if num == 2
        load 'ECGshape.mat';
        L = length(y);
        y = spline(0:1/L:(L-1)/L,y,x);
        y = y - mean(y);
        y = y/sqrt(sum(abs(y).^2)/length(y));
    else if num ==3
        load 'ECGshape2.mat';
        L = length(y);
        y = spline(0:1/L:(L-1)/L,y,x);
        y = y - mean(y);
        y = y/sqrt(sum(abs(y).^2)/length(y));
        else if num == 4
            loc = find(x<=0.4);
            y(loc) = x(loc);
            loc = find(x>0.4 & x<=0.6);
            y(loc) = 2-4*x(loc);
            loc = find(x>0.6);
            y(loc) = x(loc)-1;
            y = y - mean(y);
            y = y/sqrt(sum(abs(y).^2)/length(y));
            else
                loc = find(x<=0.45);
                y(loc) = x(loc);
                loc = find(x>0.45 & x<=0.55);
                y(loc) = 4.5-9*x(loc);
                loc = find(x>0.55);
                y(loc) = x(loc)-1;
                y = y - mean(y);
                y = y/sqrt(sum(abs(y).^2)/length(y));
            end
        end
    end
end