function CIRCLE(x_c,y_c,r)
% This function plots a circle using its center and radius.

theta = 0 : pi/50 : 2*pi;
x = r * cos(theta) + x_c;
y = r * sin(theta) + y_c;
plot(x,y);

end

