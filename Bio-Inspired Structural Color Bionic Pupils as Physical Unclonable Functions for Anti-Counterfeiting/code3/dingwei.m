% function [circleiris,circlepupil,imagewithnoise]=dingwei(eyeimage_filename,write)

%normalisation parameters
radial_res = 100;
angular_res = 240;
% with these settings a 9600 bit iris template is
% created

%feature encoding parameters
nscales=1;
minWaveLength=18;
mult=1; % not applicable if using nscales = 1
sigmaOnf=0.5;


% IO
eyeimage = imread(eyeimage_filename); 
eyeimage = rgb2gray(eyeimage);
savefile = [eyeimage_filename(1:end-4) ,'-houghpara.mat'];
[stat,mess]=fileattrib(savefile);
if stat == 1
    % if this file has been processed before
    % then load the circle parameters and
    % noise information for that file.
    load(savefile);
else
    cd '../myseg' %调用seg里面的虹膜区域分割，而且隔离噪声区域
    % if this file has not been processed before then perform automatic segmentation and save the results to a file
    % 調用segmentris得到 circleiris, circlepupil, imagewithnoise
    [circleiris,circlepupil,imagewithnoise] = segmentiris(eyeimage);
    cd '../code3'
    save(savefile,'circleiris','circlepupil','imagewithnoise');
    
end

%WRITE NOISE IMAGE

% 设置为整数 8字节
% imagewithnoise的形状是[x，y]的二维矩阵，其中有噪声的地方被标记为NaN
imagewithnoise2 = uint8(imagewithnoise);
imagewithcircles = uint8(eyeimage);

%get pixel coords for circle around iris
[x,y] = circlecoords([circleiris(2),circleiris(1)],circleiris(3),size(eyeimage));
ind2 = sub2ind(size(eyeimage),double(y),double(x)); 

%get pixel coords for circle around pupil
[xp,yp] = circlecoords([circlepupil(2),circlepupil(1)],circlepupil(3),size(eyeimage));
% subjective to index 下标转换索引， 得到一个数组，数组对应了iris/pupil的坐标位置
% 那这个时候ind(ex)1 
ind1 = sub2ind(size(eyeimage),double(yp),double(xp));

% 此处就是把这些iris和Pupil位置都设置为"白色"
% Write noise regions
imagewithnoise2(ind2) = 255;
imagewithnoise2(ind1) = 255;
% Write circles overlayed
imagewithcircles(ind2) = 255;
imagewithcircles(ind1) = 255;

% 要存图了
% w = cd;
% cd(DIAGPATH);
pos = findstr(eyeimage_filename,'\');
posdot = findstr(eyeimage_filename,'.');
l = length(pos);
addpos = pos(l);

final_segmented = [eyeimage_filename(1:addpos),'segmented-',eyeimage_filename(addpos+1:posdot),'.jpeg'];


if write   
final_noise = [eyeimage_filename(1:addpos),'noise-',eyeimage_filename(addpos+1:posdot),'.jpeg'];    
imwrite(imagewithnoise2,final_noise,'jpeg');    
%write the *-gabor_original.jpg
writeoriginal(circleiris,circlepupil,eyeimage,eyeimage_filename,nscales, minWaveLength, mult, sigmaOnf);
end 
