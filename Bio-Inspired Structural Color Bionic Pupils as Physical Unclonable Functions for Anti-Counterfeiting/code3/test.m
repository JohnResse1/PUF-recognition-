% I2=imread('image004.jpg');
 I2=imread('0023_006.bmp');

figure,imshow(I2);

eI=edge(I2,'canny', 0.2);
figure,imshow(eI);

[y0detect,x0detect,Accumulator] = houghcircle(eI,45,4);

figure,imshow(I2);
hold;
for i=1:length(y0detect)
    plot(x0detect,y0detect,'.r');
    hold on;
end

figure,imshow(Accumulator,[]);
[r,c]=size(I2);

M = circle( c,r,x0detect,y0detect,45);
figure,imshow(M,[]);

outI=M.*double(I2);
figure,imshow(outI,[]);

outI2=(1-M).*double(I2);
figure,imshow(outI2,[]);
