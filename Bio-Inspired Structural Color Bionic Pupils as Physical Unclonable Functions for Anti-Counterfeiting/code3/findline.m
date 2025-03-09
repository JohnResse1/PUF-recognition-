% findline - returns the coordinates of a line in an image using the
% linear Hough transform and Canny edge detection to create
% the edge map.（利用线性Hough变换和Canny边缘探测得到的线上各点的坐标）
%
% Usage: 
% lines = findline(image)
%
% Arguments:
%	image   - the input image
%
% Output:
%	lines   - parameters of the detected line in polar form


function lines = findline(image)

% [I2 or] = canny(image, 2, 1, 0.00, 1.00);%function [gradient, or] = canny(im, sigma, scaling, vert, horz),I2为灰度阶梯图
% 
% I3 = adjgamma(I2, 1.9);%function newim = adjgamma(im, g)
% I4 = nonmaxsup(I3, or, 1.5);%im = nonmaxsup(inimage, orient, radius)
% edgeimage = hysthresh(I4, 0.20, 0.15);%function bw = hysthresh(im, T1, T2)

% generate the edge image using toolbox
edgeimage = edge(image,'canny');

theta = (0:179)';
[R, xp] = radon(edgeimage, theta);

maxv = max(max(R));

if maxv > 25
    i = find(R == max(max(R)));
else
    lines = [];
    return;
end

[foo, ind] = sort(-R(i));
u = size(i,1);
k = i(ind(1:u));
[y,x]=ind2sub(size(R),k);
t = -theta(x)*pi/180;
r = xp(y);

lines = [cos(t) sin(t) -r];

cx = size(image,2)/2-1;
cy = size(image,1)/2-1;
lines(:,3) = lines(:,3) - lines(:,1)*cx - lines(:,2)*cy;
