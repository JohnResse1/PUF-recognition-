% shiftbits - function to shift the bit-wise iris patterns in order to provide the best match
% each shift is by two bit values and left to right, since one pixel value in the
% normalised iris pattern gives two bit values in the template
% also takes into account the number of scales
% used(像素偏移以进行归一化后的图像匹配识别，消除人眼球旋转早晨的误差？）
%
% Usage: 
% [template, mask] = createiristemplate(eyeimage_filename)
%
% Arguments:
%	template        - the template to shift
%   noshifts        - number of shifts to perform to the right, a negative
%                     value results in shifting to the left
%   nscales         - number of filters used for encoding, needed to
%                     determine how many bits to move in a shift
%
% Output:
%   templatenew     - the shifted template
%


function templatenew = shiftbits(template, noshifts,nscales)

templatenew = zeros(size(template));

width = size(template,2);
s = round(2*nscales*abs(noshifts));
p = round(width-s);

if noshifts == 0
    templatenew = template;
    
    % if noshifts is negatite then shift towards the left
elseif noshifts < 0
    
    x=1:p;
    
    templatenew(:,x) = template(:,s+x);
    
    x=(p + 1):width;
    
    templatenew(:,x) = template(:,x-p);
    
else
    
    x=(s+1):width;
    
    templatenew(:,x) = template(:,x-s);
    
    x=1:s;
    
    templatenew(:,x) = template(:,p+x);
    
end