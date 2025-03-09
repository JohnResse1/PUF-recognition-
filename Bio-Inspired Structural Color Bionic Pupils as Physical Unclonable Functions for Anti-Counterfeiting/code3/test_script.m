% 指定文件夹路径
folderPath ='dataset\train';

% 获取该文件夹及其所有子文件夹的路径
allFiles = dir(fullfile(folderPath, '**', '*'));

% 遍历所有文件
for k = 1:length(allFiles)
    % 获取文件名
    fileName = allFiles(k).name;
    
    % 排除当前和父文件夹
    if ~strcmp(fileName, '.') && ~strcmp(fileName, '..')
        fullPath = fullfile(allFiles(k).folder, fileName);
        if ~allFiles(k).isdir % 只处理文件
            % 检查文件是否以.jpg结尾
            [~, ~, ext] = fileparts(fileName);
            if ~strcmpi(ext, '.jpg')
                % 删除文件
                delete(fullPath);
            end
        end
    end
end

folderPath = 'dataset\different';
% 获取该文件夹及其所有子文件夹中的所有.jpg文件
jpgFiles = dir(fullfile(folderPath, '**', '*.jpg'));

% 遍历所有找到的.jpg文件
fileCounter = containers.Map; % 创建一个映射来记录每个子文件夹的计数器

for k = 1:length(jpgFiles)
    % 获取当前文件的完整路径
    fullPath = fullfile(jpgFiles(k).folder, jpgFiles(k).name);
    
    % 获取子文件夹名称
    [~, subFolderName] = fileparts(jpgFiles(k).folder);
    
    % 更新计数器
    if isKey(fileCounter, subFolderName)
        fileCounter(subFolderName) = fileCounter(subFolderName) + 1;
    else
        fileCounter(subFolderName) = 1;
    end
    
    % 创建新的文件名，序号格式为两位数
    newFileName = sprintf('%s_%02d.jpg', subFolderName, fileCounter(subFolderName));
    
    % 设置新的完整路径
    newFullPath = fullfile(jpgFiles(k).folder, newFileName);
    
    % 重命名文件
    try
    movefile(fullPath, newFullPath);
    end
end

% 设置要遍历的主文件夹和目标文件夹
sourceFolder = 'D:\20241003\PUF-png - 副本\dataset';
targetFolder = 'dataset\different';

% 获取主文件夹下的所有子文件夹
subFolders = dir(fullfile(sourceFolder, '**', '')); 
subFolders = subFolders([subFolders.isdir]); % 仅保留文件夹

% 初始化计数器
fileCounter = 1;

% 遍历每个子文件夹
for i = 1:length(subFolders)
    if ~strcmp(subFolders(i).name, '.') && ~strcmp(subFolders(i).name, '..')
        % 获取子文件夹的完整路径
        subFolderPath = fullfile(subFolders(i).folder, subFolders(i).name);
        
        % 查找子文件夹中的png文件
        pngFiles = dir(fullfile(subFolderPath, '*.png'));
        
        % 如果找到png文件，则复制第一张到目标文件夹并重命名
        if ~isempty(pngFiles)
            sourceFile = fullfile(subFolderPath, pngFiles(1).name);
            
            % 读取图像
            img = imread(sourceFile);
            
            % 调整图像大小为300x300
            resizedImg = imresize(img, [300 300]);
            
            % 生成新的文件名
            newFileName = sprintf('diff_%03d.jpg', fileCounter);
            targetFile = fullfile(targetFolder, newFileName);
            
            % 保存调整大小后的图像为jpg格式
            imwrite(resizedImg, targetFile);
            
            % 更新计数器
            fileCounter = fileCounter + 1;
        end
    end
end

% 设置路径和参数
aa = '06';
sourceFolder = ['D:\20241003\PUF-png - 副本\dataset\' aa '\'];
targetFolder = ['dataset\train\' aa '\'];
numImagesToGenerate = 100; % 总共生成的图片数量
% % 设置路径和参数
% aa = 'same';
% sourceFolder = 'D:\20241003\PUF-png - 副本\1';
% targetFolder = 'dataset\same';

numImagesToGenerate = 100; % 总共生成的图片数量
% 获取所有图片文件
imageFiles = dir(fullfile(sourceFolder, '*.png'));
numImages = length(imageFiles);

% 确保目标文件夹存在
if ~exist(targetFolder, 'dir')
    mkdir(targetFolder);
end

imageCounter = 0;

% 随机生成图像
for i = 1:numImagesToGenerate
    % 随机选择一张图片
    randomImageIndex = randi(numImages);
    img = imread(fullfile(sourceFolder, imageFiles(randomImageIndex).name));
    
    img = imresize(img, [300 300]);
    % 随机生成旋转角度和缩放比例
    randomScale = 0.99 + rand() * 0.02; % 在0.99到1.01之间
    randomAngle = rand() * 1; % 在0到5度之间
    
    % 旋转图片
    rotatedImg = imrotate(img, randomAngle, 'bilinear', 'crop');
    
    % 缩放图片
    scaledImg = imresize(rotatedImg, randomScale);
    
    % 如果缩放后的图片小于300x300，则需要补充边界
    [height, width, ~] = size(scaledImg);
    if height < 300 || width < 300
        % 使用 padarray 对图像进行边缘补充
        scaledImg = padarray(scaledImg, [max(0, 300 - height), max(0, 300 - width)], 'replicate', 'post');
    end
    
    % 将图片裁剪成300x300，裁剪从中心区域进行
    centerX = floor(size(scaledImg, 2) / 2);
    centerY = floor(size(scaledImg, 1) / 2);
    rect = [centerX - 150, centerY - 150, 299, 299]; % [x, y, width, height]
    croppedImg = imcrop(scaledImg, rect);
    
    % 生成新的文件名
    newFileName = sprintf('%s_%02d.jpg', aa, imageCounter);
    imwrite(croppedImg, fullfile(targetFolder, newFileName));
    
    % 更新计数器
    imageCounter = imageCounter + 1;
end

% 设置源文件夹和目标文件夹
sourceFolder = 'D:\20241003\PUF-png - 副本\dataset\新建文件夹';
targetFolder = 'dataset\different';

% 获取源文件夹中的所有png文件
pngFiles = dir(fullfile(sourceFolder, '*.png'));

% 初始化计数器
fileCounter = 1;

% 遍历每个png文件
for i = 1:length(pngFiles)
    % 获取当前png文件的完整路径
    sourceFile = fullfile(pngFiles(i).folder, pngFiles(i).name);
    
    % 读取图像
    img = imread(sourceFile);
    
    % 调整图像大小为300x300
    resizedImg = imresize(img, [300 300]);
    
    % 生成新的文件名，以序号重命名
    newFileName = sprintf('diff_%03d.jpg', fileCounter);
    targetFile = fullfile(targetFolder, newFileName);
    
    % 保存调整大小后的图像
    imwrite(resizedImg, targetFile);
    
    % 更新计数器
    fileCounter = fileCounter + 1;
end
