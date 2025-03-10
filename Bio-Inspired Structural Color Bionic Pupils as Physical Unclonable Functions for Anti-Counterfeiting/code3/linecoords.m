% linecoords - returns the x y coordinates of positions along a
% line（返回一条线的上各点的x,y坐标值）
%
% Usage: 
% [x,y] = linecoords(lines, imsize)
%
% Arguments:
%	lines       - an array containing parameters of the line in
%                 form
%   imsize      - size of the image, needed so that x y coordinates
%                 are within the image boundary
%
% Output:
%	x           - x coordinates
%	y           - corresponding y coordinates


function [x,y] = linecoords(lines, imsize)

xd = [1:imsize(2)];
yd = (-lines(3) - lines(1)*xd ) / lines(2);

coords = find(yd>imsize(1));
yd(coords) = imsize(1);
coords = find(yd<1);
yd(coords) = 1;

x = int32(xd);
y = int32(yd);   