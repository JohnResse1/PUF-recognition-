%  （对经过归一化的图像进行特征提取并编码得到E0，E0为复数矩阵）
[E0,filtersum] = gaborconvolve(polar_array, nscales, minWaveLength, mult, sigmaOnf);
% save E0 E0
f0=double(E0{1,1});
f = abs(f0);  
%%
% 编码，产生特征模板
length1 = size(polar_array,2)*2*nscales;
template = zeros(size(polar_array,1), length1);
length2 = size(polar_array,2);%列数
h = 1:size(polar_array,1);%h为行数
%create the iris template
mask = zeros(size(template));
for k=1:nscales%循环，给几个滤波器，1
    E1 = E0{k};
%     
%     % WRITE NORMALISED PATTERN, AND NOISE PATTERN
%     w = cd;
%     cd(testpath);
%     imwrite(E1,[eyeimage_filename,'-gabor_normal.jpg'],'jpg');
%     %imwrite(noise_array,[eyeimage_filename,'-polarnoise.jpg'],'jpg');
%     cd(w); 
    %Phase quantisation ： gabor过的图象被二值化
    H1 = real(E1) > 0;%实部可以对图像进行平滑滤波，虚部可以用来边缘检测
    H2 = imag(E1) > 0;
    % if amplitude is close to zero then
    % phase data is not useful, so mark off
    % in the noise mask
    H3 = abs(E1) < 0.0001;
    for i=0:(length2-1)       
        ja = double(2*nscales*(i));
        %construct the biometric template 生物特征模板
        template(h,ja+(2*k)-1) = H1(h, i+1);
        template(h,ja+(2*k)) = H2(h,i+1);
        %create noise mask 噪声掩码
        mask(h,ja+(2*k)-1) = noise_array(h, i+1) | H3(h, i+1);
        mask(h,ja+(2*k)) =   noise_array(h, i+1) | H3(h, i+1);
    end
end 