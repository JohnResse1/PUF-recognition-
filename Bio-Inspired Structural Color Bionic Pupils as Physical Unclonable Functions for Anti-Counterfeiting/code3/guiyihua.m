% perform normalisation
savefile = [eyeimage_filename(1:end-4),'-houghpara.mat'];
load(savefile); %'circleiris','circlepupil','imagewithnoise'
[polar_array,noise_array] = normaliseiris(imagewithnoise, circleiris(2),...
    circleiris(1), circleiris(3), circlepupil(2), circlepupil(1), circlepupil(3),eyeimage_filename, radial_res, angular_res,write);

%
aa=polar_array;
bb=noise_array;
% figure;subplot(121);imshow(aa);subplot(122);imshow(bb)
%

pos = findstr(eyeimage_filename,'\');
posdot = findstr(eyeimage_filename,'.');
l = length(pos);
addpos = pos(l);
if write
final_polar = [eyeimage_filename(1:addpos),'polar-',eyeimage_filename(addpos+1:posdot),'.jpeg'];
final_polarnoise = [eyeimage_filename(1:addpos),'polarnoise-',eyeimage_filename(addpos+1:posdot),'.jpeg'];
imwrite(polar_array,final_polar,'jpeg');
imwrite(noise_array,final_polarnoise,'jpeg');
% cd(w);
end
