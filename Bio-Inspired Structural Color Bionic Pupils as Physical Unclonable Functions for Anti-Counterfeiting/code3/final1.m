
function [result,time] = final1(testimage)

t0 = clock;
hmthresh = 0.3;
write = 0;
InputPath='.\0023\';
FileName=dir(strcat(InputPath,'*.bmp'));
NumFile=length(FileName);
hd=zeros(1,NumFile);
[templatetest, masktest] = createiristemplate(testimage,write);

for i=1:NumFile
    tempFileName=FileName(i).name;
    ImPath=strcat(InputPath,tempFileName);
[template, mask] = createiristemplate(ImPath,write);
hd(i) = gethammingdistance(templatetest, masktest, template, mask, 4);
   if hd(i) < hmthresh 
       result =FileName(i).name;
      break;
   end 
end

if i== NumFile
    k = find(hd==min(hd));
result = FileName(k).name;
else
    result = 'oo, No match found!';
end

time = etime(clock, t0);