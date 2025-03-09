% gaborconvolve - function for convolving each row of an image with 1D
% log-Gabor filters（利用一维盖伯滤波器对图像进行处理）
% Usage: 
% [template, mask] = createiristemplate(eyeimage_filename)
%
% Arguments:
%   im              - the image to convolve
%   nscale          - number of filters to use
%   minWaveLength   - wavelength of the basis filter
%   mult            - multiplicative factor between each filter
%   sigmaOnf        - Ratio of the standard deviation of the Gaussian describing
%                     the log Gabor filter's transfer function in the frequency domain to the filter center frequency.
% Output:
%   E0 - a 1D cell array of complex valued comvolution results
function [EO, filtersum] = gaborconvolve(im, nscale, minWaveLength, mult, sigmaOnf)%im:归一化后的polar-array

[rows,cols] = size(im);		
%filtersum = zeros(1,size(im,2));

EO = cell(1, nscale);          % Pre-allocate cell array

ndata = cols;
if mod(ndata,2) == 1             % If there is an odd No of data points 
    ndata = ndata-1;               % throw away the last one.
end

logGabor  = zeros(1,ndata);
result = zeros(rows,ndata);
filtersum = logGabor;

radius =  [0:fix(ndata/2)]/fix(ndata/2)/2;  % Frequency values 0 - 0.5
radius(1) = 1;
wavelength = minWaveLength;        % Initialize filter wavelength.
for s = 1:nscale                  % For each scale.  
    
    % Construct the filter - first calculate the radial filter component.
    fo = 1.0/wavelength;                  % Centre frequency of filter.
    rfo = fo/0.5;                         % Normalised radius from centre of frequency plane 
    % corresponding to fo.
    logGabor(1:ndata/2+1) = exp((-(log(radius/fo)).^2) / (2 * log(sigmaOnf)^2));  
    logGabor(1) = 0;  
    filter = logGabor;
    filtersum = filtersum+filter;
    % for each row of the input image, do the convolution, back transform
    for r = 1:rows	% For each row
        signal = im(r,1:ndata);%每一行所有列给到signal
        imagefft = fft( signal );%经过fft变到频域
        result(r,:) = ifft(imagefft .* filter);%点乘滤波器后变回去得到滤波结果
    end
    % save the ouput for each scale
    EO{s} = result;%滤波结果放到E0
    wavelength = wavelength * mult;       % Finally calculate Wavelength of next filter
end                                     % ... and process the next scale
filtersum = fftshift(filtersum);