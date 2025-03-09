%  ���Ծ�����һ����ͼ�����������ȡ������õ�E0��E0Ϊ��������
[E0,filtersum] = gaborconvolve(polar_array, nscales, minWaveLength, mult, sigmaOnf);
% save E0 E0
f0=double(E0{1,1});
f = abs(f0);  
%%
% ���룬��������ģ��
length1 = size(polar_array,2)*2*nscales;
template = zeros(size(polar_array,1), length1);
length2 = size(polar_array,2);%����
h = 1:size(polar_array,1);%hΪ����
%create the iris template
mask = zeros(size(template));
for k=1:nscales%ѭ�����������˲�����1
    E1 = E0{k};
%     
%     % WRITE NORMALISED PATTERN, AND NOISE PATTERN
%     w = cd;
%     cd(testpath);
%     imwrite(E1,[eyeimage_filename,'-gabor_normal.jpg'],'jpg');
%     %imwrite(noise_array,[eyeimage_filename,'-polarnoise.jpg'],'jpg');
%     cd(w); 
    %Phase quantisation �� gabor����ͼ�󱻶�ֵ��
    H1 = real(E1) > 0;%ʵ�����Զ�ͼ�����ƽ���˲����鲿����������Ե���
    H2 = imag(E1) > 0;
    % if amplitude is close to zero then
    % phase data is not useful, so mark off
    % in the noise mask
    H3 = abs(E1) < 0.0001;
    for i=0:(length2-1)       
        ja = double(2*nscales*(i));
        %construct the biometric template ��������ģ��
        template(h,ja+(2*k)-1) = H1(h, i+1);
        template(h,ja+(2*k)) = H2(h,i+1);
        %create noise mask ��������
        mask(h,ja+(2*k)-1) = noise_array(h, i+1) | H3(h, i+1);
        mask(h,ja+(2*k)) =   noise_array(h, i+1) | H3(h, i+1);
    end
end 