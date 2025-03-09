% ADJGAMMA - Adjusts image gamma.（调整图像的伽马值）
%
% function g = adjgamma(im, g)


function newim = adjgamma(im, g)

    if g <= 0
	error('Gamma value must be > 0');
    end

    if isa(im,'uint8');
	newim = double(im);
    else 
	newim = im;
    end
    	
    % rescale range 0-1
    newim = newim-min(min(newim));
    newim = newim./max(max(newim));
    
    newim =  newim.^(1/g);   % Apply gamma function


