function [M] = circle( x,y,cx,cy,r)

[x,y]= meshgrid(-(cx-1):(x-cx),-(cy-1):(y-cy));
surf(x,y,x.^2+y.^2);
M=((x.^2+y.^2)<=r.^2);

