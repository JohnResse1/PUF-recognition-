function writeoriginal(circleiris,circlepupil,eyeimage,eyeimage_filename,nscales, minWaveLength, mult, sigmaOnf)
   global DIAGPATH
    eyeimage1 = eyeimage;
    circleirisd = double(circleiris);
    circlepupild = double(circlepupil);
    
    row_iris = circleirisd(1);
    clo_iris = circleirisd(2);
    rad_iris = circleirisd(3);
    l_iris   = clo_iris-rad_iris;
    r_iris   = clo_iris+rad_iris;
    t_iris   = row_iris-rad_iris;
    b_iris   = row_iris+rad_iris;
    
    row_pupil= circlepupild(1);
    clo_pupil= circlepupild(2);
    rad_pupil= circlepupild(3);
    l_pupil   = clo_pupil-rad_pupil;
    r_pupil   = clo_pupil+rad_pupil;
    t_pupil   = row_pupil-rad_pupil;
    b_pupil   = row_pupil+rad_pupil;
    
    %大于边界的位置都置为0，使得留下的位置是我们感兴趣的区域，而为0的位置则不处理。也就是图像的背底。
    %original = zeros(t_iris:b_iris,l_iris:r_iris);
    [row,clo]=size(eyeimage);
    for i=1: row
        for k=1: clo
            dis1 = ((i-row_iris)^2+(k-clo_iris)^2)^0.5;
            dis2 = ((i-row_pupil)^2+(k-clo_pupil)^2)^0.5;
            % 对于eyeimage考虑大于iris和小于pupil之间的部分，置为0
            if (dis1 > rad_iris) || (dis2 < rad_pupil)
                eyeimage(i,k) = 0;
            end
            % 对于eyeimage1考虑大于iris的部分
            if (dis1 > rad_iris)
                eyeimage1(i,k) = 0;
            end            
        end
    end

    
    original = eyeimage(t_iris:b_iris,l_iris:r_iris);
    original1 = eyeimage1(t_iris:b_iris,l_iris:r_iris);   
      

   % Garbor滤波边缘检测
   % 对于刚刚我们置0的明显边缘图像提取边缘位置
   for k=1:nscales  
       [E0 filtersum0] = gaborconvolve(double(original), nscales, minWaveLength, mult, sigmaOnf);
       [E1 filtersum1] = gaborconvolve(double(original1), nscales, minWaveLength, mult, sigmaOnf);
       

       % max是有的就留下，没有的还是0
       writeimage0=abs(real(E0{k}));
         Max0=max(max(writeimage0));
         
             [r,c]=size(writeimage0);
     mask = ones(r,c);

    for j=1:r
        for g = 1:c
            dis = ((j+t_iris-row_iris)^2+(g+l_iris-clo_iris)^2)^0.5;
            if dis > rad_iris
                mask(j,g)= 0;
            end
        end
    end
         % 这里提取了iris轮廓
         writeimage0=adjgamma(Max0-writeimage0,0.01).*mask;
       % 同理这里提取了pupils
       writeimage1=(real(E1{k}));


pos = findstr(eyeimage_filename,'\');
posdot = findstr(eyeimage_filename,'.');
l = length(pos);
addpos = pos(l);
final_gabor_original = [eyeimage_filename(1:addpos),'gabor_original-',eyeimage_filename(addpos+1:posdot),'.jpeg'];     
final_gabor_lins = [eyeimage_filename(1:addpos),'gabor_lins-',eyeimage_filename(addpos+1:posdot),'.jpeg']; 

imwrite(writeimage0,final_gabor_original,'jpeg');
imwrite(writeimage1,final_gabor_lins,'jpeg');


%  
%     testpath = [DIAGPATH , '\testDIA'];
%     w = cd;
%     cd(testpath); 

%     imwrite(writeimage1,[eyeimage_filename,'-gabor_original1.jpg'],'jpg');
%     cd(w);
    
    
    
   end