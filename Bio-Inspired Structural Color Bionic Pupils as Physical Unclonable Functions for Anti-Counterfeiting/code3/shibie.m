
% [templatetest, masktest] = createiristemplate(testimage,write);
t0 = clock;
FileName=dir(strcat(InputPath,'*.bmp'));
NumFile=length(FileName);
hd=zeros(1,NumFile);
write=0;
txtname=['HanmingDist_',filenamestr,'.txt'];
if ~exist(txtname)
    file_id=fopen(txtname,'a+');
end
fid=fopen(txtname,'w');
for i=1:NumFile
    tempFileName=FileName(i).name;
    ImPath=strcat(InputPath,tempFileName);
[template, mask] = createiristemplate(ImPath,write);
hd(i) = gethammingdistance(templatetest, masktest, template, mask, 4);
    fprintf(fid,'%s\t%f\n',FileName(i).name,hd(i));
%    if hd(i) < hmthresh 
%        result =FileName(i).name;
%       break;
%    end 
end
fclose(fid);
%找最小距离的作为识别结果
[~,k]=min(hd);k=k(1);
result = FileName(k).name;
time = etime(clock, t0);