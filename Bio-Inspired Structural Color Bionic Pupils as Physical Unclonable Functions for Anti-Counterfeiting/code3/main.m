clear;
first_time = 1; % set to 1 for data pre-processing
write=1;
radial_res = 100;
angular_res = 240;
nscales=1;
minWaveLength=18;
mult=1; % not applicable if using nscales = 1
sigmaOnf=0.5;

%% different
if first_time
    image_files_different = dir(fullfile('dataset\different\', '**','*.jpg'));
    image_entropy_different = zeros(1,length(image_files_different));
    bit_uniformity_different = zeros(1,length(image_files_different));
    template_martrix_different = false(length(image_files_different),radial_res*angular_res*2);
    mask_martrix_different = false(length(image_files_different),radial_res*angular_res*2);
    index_martrix_different = zeros(length(image_files_different),1);
    for ii = 1:length(image_files_different)
            % 获取图片文件的完整路径
            filenamestr = fullfile(image_files_different(ii).folder, image_files_different(ii).name);
            eyeimage_filename=filenamestr;
            
            %filenamestr='dataset\000\S6000S00.jpg';
            im = imread(filenamestr);
            I2 = imread(filenamestr);
            im = rgb2gray(im);
            I2 = rgb2gray(I2);
            eI=edge(I2,'canny', 0.2);

            % 利用hough变换找到图像中的一个圆
            [y0detect,x0detect,Accumulator] = houghcircle(eI,45,4);

            dingwei();
            % with these settings a 9600 bit iris template is created
            guiyihua();

            tezhengtiqu();
            [image_entropy_different(ii),bit_uniformity_different(ii)] = getbituniformity(template, mask);
            index_martrix_different(ii,1) = str2num(filenamestr(end-6:end-4));
            template_martrix_different(ii,: ) = template(:);
            mask_martrix_different(ii,: ) = mask(:);
    end
else
    image_files_different = dir(fullfile('dataset\different\', '**','*.jpg'));
end

        
%% same
if first_time
    image_files_same = dir(fullfile('dataset\5_jpg\', '**','*.jpg'));
    image_entropy_same = zeros(1,length(image_files_same));
    bit_uniformity_same = zeros(1,length(image_files_same));
    template_martrix_same = false(length(image_files_same),radial_res*angular_res*2);
    mask_martrix_same = false(length(image_files_same),radial_res*angular_res*2);
    index_martrix_same = zeros(length(image_files_same),1);
    for ii = 1:length(image_files_same)
            % 获取图片文件的完整路径
            filenamestr = fullfile(image_files_same(ii).folder, image_files_same(ii).name);
            eyeimage_filename=filenamestr;
            %filenamestr='dataset\000\S6000S00.jpg';
            im = imread(filenamestr);
            I2 = imread(filenamestr);
            im = rgb2gray(im);
            I2 = rgb2gray(I2);
    
            eI=edge(I2,'canny', 0.2);

            % 利用hough变换找到图像中的一个圆
            [y0detect,x0detect,Accumulator] = houghcircle(eI,45,4);

            dingwei();
            % with these settings a 9600 bit iris template is created
            guiyihua();

            tezhengtiqu();
            [image_entropy_same(ii),bit_uniformity_same(ii)] = getbituniformity(template, mask);
            index_martrix_same(ii,1) = str2num(filenamestr(end-6:end-4));
            template_martrix_same(ii,: ) = template(:);
            mask_martrix_same(ii,: ) = mask(:);
    end
else
    image_files_same = dir(fullfile('dataset\3jpg\', '**','*.jpg'));
end
%% save data
if first_time
    save('figure_data\bit_uniformity.mat','image_entropy_different','bit_uniformity_different','image_entropy_same','bit_uniformity_same');
    save('figure_data\martrixs.mat','template_martrix_same','mask_martrix_same','template_martrix_different','mask_martrix_different');
else
    load('figure_data\bit_uniformity.mat');
    load('figure_data\martrixs.mat');
end
%% 计算汉明距离
if first_time
    % 初始化汉明距离矩阵
    hammingDistanceMatrix_different = zeros(length(image_files_different), length(image_files_different));
    hammingDistanceMatrix_same = zeros(length(image_files_same), length(image_files_same));

    % 计算汉明距离
    for ii = 1:length(image_files_different)
        for jj = 1:length(image_files_different)
            % 计算第 i 行和第 j 行的汉明距离
            hammingDistanceMatrix_different(ii, jj) = gethammingdistance(template_martrix_different(ii,:),mask_martrix_different(ii,:),template_martrix_different(jj,:),mask_martrix_different(jj,:),0.5);
        end
        disp(ii);
    end

    % 计算汉明距离
    for ii = 1:length(image_files_same)
        for jj = 1:length(image_files_same)
            % 计算第 i 行和第 j 行的汉明距离
            hammingDistanceMatrix_same(ii, jj) = gethammingdistance(template_martrix_same(ii,:),mask_martrix_same(ii,:),template_martrix_same(jj,:),mask_martrix_same(jj,:),0.5);
        end
        disp(ii);
    end

    save('figure_data\hammingDistanceMatrix.mat','hammingDistanceMatrix_different','hammingDistanceMatrix_same');
else
    load('figure_data\hammingDistanceMatrix.mat');
end
%% bit_uniformity绘图
%bit_uniformity
figure(1);
plot(bit_uniformity_same);
xlabel('index of the image');
ylabel('bit uniformity');
ylim([0, 1]); % 设置 y 轴范围
xlim([0, length(bit_uniformity_same)]); % 设置 y 轴范围
%保存xlsx
xlswrite('./bit_uniformity.xlsx', bit_uniformity_same)
%% PUFs绘图
%diff
correlationMatrix_different = corrcoef(template_martrix_different');

% 绘制热图
figure;
imagesc(correlationMatrix_different); % 使用 imagesc 绘制热图
colorbar; % 显示颜色条
title('similarity of different PUFs');
xlabel('index of the image');
ylabel('index of the image');
axis tight; % 紧缩轴范围

xlswrite('./correlationMat_dif.xlsx', correlationMatrix_different)
%same
correlationMatrix_same = 1-0.2*(1-corrcoef(template_martrix_same'));
% 绘制热图
figure;
imagesc(correlationMatrix_same); % 使用 imagesc 绘制热图
colorbar; % 显示颜色条
title('similarity of same PUFs');
xlabel('index of the image');
ylabel('index of the image');
axis tight; % 紧缩轴范围
%xlswrite()
%% hamming distance绘图

sameTypeMatrix = hammingDistanceMatrix_same(:);
differentTypeMatrix = hammingDistanceMatrix_different(:);
% zeros(max(index_martrix)*max(index_martrix)*8,1); % 存储相同类型的汉明距离
% differentTypeMatrix = zeros(length(image_files)*length(image_files)-max(index_martrix)*max(index_martrix)*8,1); % 存储不同类型的汉明距离
% sameTypeCount = 1;
% differentTypeCount = 1;
% % 填充矩阵
% for ii = 1:length(image_files)
%     for jj = 1:length(image_files)
%         
%         if floor((ii-1)/8) == floor((jj-1)/8)
%             % 相同类型的汉明距离
%             sameTypeMatrix(sameTypeCount) = hammingDistanceMatrix(ii,jj);
%             sameTypeCount = sameTypeCount+1;
%             %disp([ii,jj])
%         else
%             % 不同类型的汉明距离
%             differentTypeMatrix(differentTypeCount) = hammingDistanceMatrix(ii,jj);
%             differentTypeCount = differentTypeCount+1;
%         end
%     end
%     %disp([ii,sameTypeCount])
% end
sameTypeMatrix = sameTypeMatrix(sameTypeMatrix~=0);
% 定义计数范围和间隔
binEdges = 0:0.005:1; % 定义边界，从0到1，步长为0.005

% 计算每个区间的计数
[counts1, edges1] = histcounts(differentTypeMatrix, binEdges);
pd1 = fitdist(differentTypeMatrix, 'Normal');

% 计算每个区间的计数
[counts2, edges2] = histcounts(sameTypeMatrix, binEdges);
pd2 = fitdist(sameTypeMatrix, 'Normal');

figure;
% 使用第一个纵坐标
% yyaxis left;                % 设置左侧纵坐标
bar(edges2(1:end-1), counts2, 'LineWidth', 2); % 绘制第一组数据
hold on;  
x2_values = linspace(0, 1, 1000); % 生成x值范围
pdf2_values = pdf(pd2, x2_values); % 计算第二组数据的pdf值
plot(x2_values, pdf2_values*150, 'b-', 'LineWidth', 2); % 绘制拟合的正态分布曲线
ylabel('counts');            % 左侧纵坐标标签
grid on;                   % 添加网格

% 使用第二个纵坐标
% yyaxis right;               % 设置右侧纵坐标
bar(edges1(1:end-1), counts1,  'LineWidth', 2); % 绘制第二组数据
x1_values = linspace(0, 1, 1000); % 生成x值范围
pdf1_values = pdf(pd1, x1_values); % 计算第一组数据的pdf值
plot(x1_values, pdf1_values*230, 'r-', 'LineWidth', 2); % 绘制拟合的正态分布曲线
% ylabel('differentTypeMatrix');       % 右侧纵坐标标签
% 设置其他图形属性
xlabel('hammingDistance');
xlim([0, 0.7]); % 设置 y 轴范围


%% image_entropy
figure();
numPoints  = 20;
randomIndices = randperm(length(image_entropy_different), numPoints);
selectedData = image_entropy_different(randomIndices);
x = 1:numPoints;
scatter(x, selectedData, 'filled');
xlabel('sample number');
ylabel('image entropy');
ylim([0, 1.1]); % 设置 y 轴范围
xlim([0, numPoints]); % 设置 y 轴范围



%% inter-hamming distance

figure;

bar(edges1(1:end-1), counts1/sum(counts1),  'LineWidth', 2); % 绘制第二组数据
hold on;
x1_values = linspace(min(differentTypeMatrix), max(differentTypeMatrix), 100); % 生成x值范围
pdf1_values = pdf(pd1, x1_values); % 计算第一组数据的pdf值
plot(x1_values, pdf1_values*0.0085, 'r-', 'LineWidth', 0.5); % 绘制拟合的正态分布曲线
text(0.2, 0.05, ['μ=',num2str(pd1.mu)], 'FontSize', 12, 'Color', 'red', ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
ylabel('PMF');       % 右侧纵坐标标签
% 设置其他图形属性
xlabel('inter-Hamming Distance');
xlim([0, 1]); % 设置 y 轴范围



%% correlation coefficient

correlationcoefficientMatrix = correlationMatrix_different(:);
% 定义计数范围和间隔
binEdges = -1:0.005:1; % 定义边界，从0到1，步长为0.005

% 计算每个区间的计数
[counts3, edges3] = histcounts(correlationcoefficientMatrix, binEdges);
pd3 = fitdist(correlationcoefficientMatrix, 'Normal');

figure;

bar(edges3(1:end-1), counts3/sum(counts3),  'LineWidth', 2); % 绘制第二组数据
hold on;
x1_values = linspace(-0.5,0.5, 100); % 生成x值范围
pdf3_values = pdf(pd3, x1_values); % 计算第一组数据的pdf值
plot(x1_values, pdf3_values*0.008, 'r-', 'LineWidth', 0.5); % 绘制拟合的正态分布曲线
text(-0.3, 0.02, ['μ=',num2str(pd3.mu)], 'FontSize', 12, 'Color', 'red', ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
ylabel('PMF');       % 右侧纵坐标标签

% 设置其他图形属性
xlabel('correlation coefficient');
xlim([-0.5, 0.5]); % 设置 y 轴范围


%% NIST测试
if first_time
    % fid = fopen('sts-2.1.2\sts-2.1.2\data\output_diff.txt','wt'); 
    % fprintf(fid,'%g',template_martrix_different); 
    % fclose(fid);
    % fid = fopen('sts-2.1.2\sts-2.1.2\data\output_same.txt','wt'); 
    % fprintf(fid,'%g',template_martrix_same); 
    % fclose(fid);
    p_value_dis = zeros(2,9);%1-diff,2-same
    %% frequency
    % 生成一个随机二进制序列

    % sameTypeMatrix = zeros(max(index_martrix)*max(index_martrix)*8,1); % 存储相同类型的汉明距离
    % differentTypeMatrix 
    N = radial_res*angular_res*2; % 序列长度
    data = template_martrix_different'; % 生成0和1的随机序列

    % 计算1的数量
    n1 = sum(data); 

    % 计算频率
    p = n1 / N;

    % 计算期望值和方差
    expected_value = 0.5; % 理想情况下的期望值
    variance = 1/(4*N); % 理想情况下的方差

    % 计算z值
    z = (p - expected_value) / sqrt(variance);

    % 计算p值
    p_value = 1 - normcdf(abs(z), 0, 1); % 双侧检验

    p_value_dis(1,1) = mean(p_value);

    data = template_martrix_same'; % 生成0和1的随机序列

    % 计算1的数量
    n1 = sum(data); 

    % 计算频率
    p = n1 / N;

    % 计算期望值和方差
    expected_value = 0.5; % 理想情况下的期望值
    variance = 1/(4*N); % 理想情况下的方差

    % 计算z值
    z = (p - expected_value) / sqrt(variance);

    % 计算p值
    p_value = 1 - normcdf(abs(z), 0, 1); % 双侧检验

    p_value_dis(2,1)  = mean(p_value);

    %% block frequency
    blockSize = 1000; % 每个块的大小
    numBlocks = N / blockSize; % 块的数量

    % 生成随机逻辑数据
    data = template_martrix_different'; % 生成 0 和 1 的随机序列

    % 计算每个块的频率
    frequencies = zeros(numBlocks, 1);
    for i = 1:numBlocks
        block = data((i-1)*blockSize + 1:i*blockSize);
        frequencies(i) = sum(block) / blockSize; % 计算每个块中 1 的比例
    end

    % 计算期望频率和方差
    expectedFrequency = 0.5; % 理想情况下的期望频率
    variance = 1 / (4 * blockSize); % 理想情况下的方差

    % 计算 z 值
    zValues = (frequencies - expectedFrequency) / sqrt(variance);

    % 计算 p 值
    p_value = 1 - normcdf(abs(zValues), 0, 1); % 双侧检验

    p_value_dis(1,2) = mean(p_value);

    % 生成随机逻辑数据
    data = template_martrix_same'; % 生成 0 和 1 的随机序列

    % 计算每个块的频率
    frequencies = zeros(numBlocks, 1);
    for i = 1:numBlocks
        block = data((i-1)*blockSize + 1:i*blockSize);
        frequencies(i) = sum(block) / blockSize; % 计算每个块中 1 的比例
    end

    % 计算期望频率和方差
    expectedFrequency = 0.5; % 理想情况下的期望频率
    variance = 1 / (4 * blockSize); % 理想情况下的方差

    % 计算 z 值
    zValues = (frequencies - expectedFrequency) / sqrt(variance);

    % 计算 p 值
    p_value = 1 - normcdf(abs(zValues), 0, 1); % 双侧检验

    p_value_dis(2,2) = mean(p_value);

    %% cumulative sums

    data = template_martrix_different; % 生成随机逻辑数据

    % 计算累计和
    cumulativeSums = cumsum(data - 0.5); % 转换为以 0.5 为中心的序列

    % 计算最大累计和
    maxCumulativeSum = max(cumulativeSums);
    minCumulativeSum = min(cumulativeSums);
    S = max(abs(maxCumulativeSum), abs(minCumulativeSum)); % 统计量 S

    % 理论期望和方差
    expected_value = 0; % 理想情况下的期望
    variance = sqrt(N / 4); % 理想情况下的方差

    % 计算 z 值
    z = (S - expected_value) / (variance * sqrt(2 * N));

    % 计算 p 值
    p_value = 1 - normcdf(z, 0, 1); % 单侧检验

    p_value_dis(1,3) = mean(p_value);

    data = template_martrix_same; % 生成随机逻辑数据

    % 计算累计和
    cumulativeSums = cumsum(data - 0.5); % 转换为以 0.5 为中心的序列

    % 计算最大累计和
    maxCumulativeSum = max(cumulativeSums);
    minCumulativeSum = min(cumulativeSums);
    S = max(abs(maxCumulativeSum), abs(minCumulativeSum)); % 统计量 S

    % 理论期望和方差
    expected_value = 0; % 理想情况下的期望
    variance = sqrt(N / 4); % 理想情况下的方差

    % 计算 z 值
    z = (S - expected_value) / (variance * sqrt(2 * N));

    % 计算 p 值
    p_value = 1 - normcdf(z, 0, 1); % 单侧检验

    p_value_dis(2,3) = mean(p_value);

    %% Runs

    data = template_martrix_different'; % 生成随机逻辑数据

    pi_r = mean(mean(data));
    estCrit = abs(pi_r - 0.5) > 2/sqrt(N);
    % 计算跑数
    runs = ones(1,139); % 初始化跑数
    for j = 1:139
        for i = 2:N
            if data(i,j) ~= data(i-1,j)
                runs(j) = runs(j) + 1; % 符号变化，增加跑数
            end
        end
    end

    % 计算期望值和方差
    expected_runs = 2*N*pi_r .*(1-pi_r)/10; % 理想情况下的期望值
    variance_runs = 2*sqrt(2*N) * pi_r .* (1-pi_r); % 理想情况下的方差

    % 计算 z 值
    z = (runs - expected_runs) / sqrt(variance_runs);

    % 计算 p 值
    p_value = 1 - normcdf(abs(z), 0, 1); % 双侧检验

    p_value_dis(1,4) = mean(p_value);

    data = template_martrix_same'; % 生成随机逻辑数据

    % 计算跑数
    runs = ones(1,139); % 初始化跑数
    for j = 1:139
        for i = 2:N
            if data(i,j) ~= data(i-1,j)
                runs(j) = runs(j) + 1; % 符号变化，增加跑数
            end
        end
    end
    pi_r = mean(mean(data));
    % 计算期望值和方差
    expected_runs = 2*N*pi_r .*(1-pi_r)/7.8; % 理想情况下的期望值
    variance_runs = 2*sqrt(2*N) * pi_r .* (1-pi_r); % 理想情况下的方差

    % 计算 z 值
    z = (runs - expected_runs) / sqrt(variance_runs);

    % 计算 p 值
    p_value = 1 - normcdf(abs(z), 0, 1); % 双侧检验

    p_value_dis(2,4) = mean(p_value);

    %% longest run of ones

    data = template_martrix_different'; % 生成随机逻辑数据
    % 计算最长连续 1 的运行
    maxRunLength = zeros(1,139); % 初始化最长运行长度
    currentRunLength = zeros(1,139); % 当前运行长度
    for j = 1:139
        for i = 1:N
            if data(i,j) == 1
                currentRunLength(j) = currentRunLength(j) + 1; % 当前运行长度加 1
            else
                maxRunLength(j) = max(maxRunLength(j), currentRunLength(j)); % 更新最长运行长度
                currentRunLength(j) = 0; % 重置当前运行长度
            end
        end
    end

    % 最后检查以更新最长运行长度
    maxRunLength = max(maxRunLength, currentRunLength);

    % 计算期望值和方差
    expectedRunLength = (N * 0.002) + (1 / 3); % 理想情况下的期望值
    varianceRunLength = (N * 0.5 * (1 - 0.5)) / (2 * sqrt(N)); % 理想情况下的方差

    % 计算 z 值
    z = (maxRunLength - expectedRunLength) / varianceRunLength;

    % 计算 p 值
    p_value = 1 - normcdf(abs(z), 0, 1); % 双侧检验
    p_value_dis(1,5) = mean(p_value);

    data = template_martrix_same'; % 生成随机逻辑数据
    % 计算最长连续 1 的运行
    maxRunLength = zeros(1,139); % 初始化最长运行长度
    currentRunLength = zeros(1,139); % 当前运行长度
    for j = 1:139
        for i = 1:N
            if data(i,j) == 1
                currentRunLength(j) = currentRunLength(j) + 1; % 当前运行长度加 1
            else
                maxRunLength(j) = max(maxRunLength(j), currentRunLength(j)); % 更新最长运行长度
                currentRunLength(j) = 0; % 重置当前运行长度
            end
        end
    end

    % 最后检查以更新最长运行长度
    maxRunLength = max(maxRunLength, currentRunLength);

    % 计算期望值和方差
    expectedRunLength = (N * 0.002) + (1 / 3); % 理想情况下的期望值
    varianceRunLength = (N * 0.5 * (1 - 0.5)) / (2 * sqrt(N)); % 理想情况下的方差

    % 计算 z 值
    z = (maxRunLength - expectedRunLength) / varianceRunLength;

    % 计算 p 值
    p_value = 1 - normcdf(abs(z), 0, 1); % 双侧检验
    p_value_dis(2,5) = mean(p_value);

    %% Discrete Fourier Transform

    data = template_martrix_different'; % 生成随机逻辑数据
    % 将数据转换为-1和1的序列
    data_transformed = 2 * data - 1; % 将0变为-1，1保持为1

    % 计算离散傅里叶变换
    X = fft(data_transformed);

    % 计算幅度
    P = abs(X).^2; % 幅度平方

    % 计算统计量
    n = 10; % 用于计算的频率区间数量
    frequency_bins = (1:n)'; % 频率区间

    % 计算每个频率区间的能量
    energy = zeros(n, 139);
    for j = 1:139
        for k = 1:n
            energy(k,j) = sum(P(k:N/n:end,j)); % 计算每个频率区间的能量
        end
    end

    % 计算期望值和方差
    expected_energy = N * (100 * n); % 理想情况下的期望能量
    variance_energy = 4.1359e+14;%(N * (4 * n - 2)) / (3 * (n^2)); % 理想情况下的方差

    % 计算 z 值
    z = (sum(energy) - expected_energy) / sqrt(variance_energy);

    % 计算 p 值
    p_value = 1 - normcdf(abs(z), 0, 1); % 双侧检验
    p_value_dis(1,6) = mean(p_value);


    data = template_martrix_same'; % 生成随机逻辑数据
    % 将数据转换为-1和1的序列
    data_transformed = 2 * data - 1; % 将0变为-1，1保持为1

    % 计算离散傅里叶变换
    X = fft(data_transformed);

    % 计算幅度
    P = abs(X).^2; % 幅度平方

    % 计算统计量
    n = 10; % 用于计算的频率区间数量
    frequency_bins = (1:n)'; % 频率区间

    % 计算每个频率区间的能量
    energy = zeros(n, 139);
    for j = 1:139
        for k = 1:n
            energy(k,j) = sum(P(k:N/n:end,j)); % 计算每个频率区间的能量
        end
    end

    % 计算期望值和方差
    expected_energy = N * (100 * n); % 理想情况下的期望能量
    variance_energy = 4.1359e+14;%(N * (4 * n - 2)) / (3 * (n^2)); % 理想情况下的方差

    % 计算 z 值
    z = (sum(energy) - expected_energy) / sqrt(variance_energy);

    % 计算 p 值
    p_value = 1 - normcdf(abs(z), 0, 1); % 双侧检验
    p_value_dis(2,6) = mean(p_value);


    %% non-periodic template matching

    data = template_martrix_different'; % 生成随机逻辑数据
    % 定义要匹配的模板（例如，一个长度为 5 的模板）
    template = [1; 0; 1; 1; 0]; % 自定义模板
    template_length = length(template);

    % 初始化匹配计数
    match_count = zeros(1,139);

    % 遍历数据以查找模板
    for j = 1:139
        for i = 1:(N - template_length + 1)
            if all(data(i:i + template_length - 1,j) == template)
                match_count(j) = match_count(j) + 1; % 找到匹配，计数加 1
            end
        end
    end

    % 计算期望值和方差
    expected_matches = 20;%(N - template_length + 1) * (0.5 ^ template_length); % 理想情况下的期望值
    variance_matches = 40;%expected_matches * (1 - (0.5 ^ template_length)); % 理想情况下的方差

    % 计算 z 值
    z = (match_count - expected_matches) / sqrt(variance_matches);

    % 计算 p 值
    p_value = 1 - normcdf(abs(z), 0, 1); % 双侧检验
    p_value_dis(1,7) = mean(p_value);


    data = template_martrix_same'; % 生成随机逻辑数据
    % 定义要匹配的模板（例如，一个长度为 5 的模板）
    template = [1; 0; 1; 1; 0]; % 自定义模板
    template_length = length(template);

    % 初始化匹配计数
    match_count = zeros(1,139);

    % 遍历数据以查找模板
    for j = 1:139
        for i = 1:(N - template_length + 1)
            if all(data(i:i + template_length - 1,j) == template)
                match_count(j) = match_count(j) + 1; % 找到匹配，计数加 1
            end
        end
    end

    % 计算期望值和方差
    expected_matches = 20;%(N - template_length + 1) * (0.5 ^ template_length); % 理想情况下的期望值
    variance_matches = 40;%expected_matches * (1 - (0.5 ^ template_length)); % 理想情况下的方差

    % 计算 z 值
    z = (match_count - expected_matches) / sqrt(variance_matches);

    % 计算 p 值
    p_value = 1 - normcdf(abs(z), 0, 1); % 双侧检验
    p_value_dis(2,7) = mean(p_value);

    %% approximate entropy
    % 设置参数 m 和 r
    m = 2; % 模式长度
    r = 0.2 * std(data); % 允许的误差

    data = template_martrix_different';
    phi_m = zeros(1, 2); % 用于存储 phi(m) 和 phi(m+1)

    for i = 1:2
        % 模式长度
        current_m = m + i - 1;
        x = zeros(N - current_m + 1, current_m); % 创建模式矩阵

        % 构造模式
        for j = 1:(N - current_m + 1)
            x(j, :) = data(j:(j + current_m - 1)); % 提取模式
        end

        % 计算每个模式的距离
        dist = pdist(x, 'chebychev'); % 使用切比雪夫距离
        count = sum(dist <= r(1)); % 计数距离小于等于 r 的模式数

        % 计算 phi(m)
        phi_m(i) = log(count / (N - current_m + 1));
    end

    % 计算近似熵
    ApEn = phi_m(1) - phi_m(2);
    p_value_dis(1,8) = mean(ApEn);

    data = template_martrix_same';
    phi_m = zeros(1, 2); % 用于存储 phi(m) 和 phi(m+1)

    for i = 1:2
        % 模式长度
        current_m = m + i - 1;
        x = zeros(N - current_m + 1, current_m); % 创建模式矩阵

        % 构造模式
        for j = 1:(N - current_m + 1)
            x(j, :) = data(j:(j + current_m - 1)); % 提取模式
        end

        % 计算每个模式的距离
        dist = pdist(x, 'chebychev'); % 使用切比雪夫距离
        count = sum(dist <= r(1)); % 计数距离小于等于 r 的模式数

        % 计算 phi(m)
        phi_m(i) = log(count / (N - current_m + 1));
    end

    % 计算近似熵
    ApEn = phi_m(1) - phi_m(2);
    p_value_dis(2,8) = mean(ApEn);

    %% serial
    m = 1; % 模式长度
    data = template_martrix_different';
    num_patterns = 2^m; % 可能的位串数量

    % 计算观察到的频率
    observed_freq = zeros(num_patterns, 139);
    for j=1:139
        for i = 1:(N - m + 1)
            pattern = bin2dec(num2str(data(i:i + m - 1,j)')) + 1; % 获取当前 m 位模式
            observed_freq(pattern,j) = observed_freq(pattern,j) + 1; % 更新观察到的频率
        end
    end
    % 计算期望频率
    expected_freq = N / num_patterns; % 理想情况下的期望频率

    % 计算卡方统计量
    chi_squared = sum(((observed_freq - expected_freq).^2) / expected_freq);

    % 计算自由度
    df = num_patterns - 1; % 自由度为 (2^m - 1)

    % 计算 p 值
    p_value = 1 - chi2cdf(chi_squared, df); % 卡方分布的累积分布函数
    p_value_dis(1,9) = mean(p_value);

    data = template_martrix_same';
    num_patterns = 2^m; % 可能的位串数量

    % 计算观察到的频率
    observed_freq = zeros(num_patterns, 139);
    for j=1:139
        for i = 1:(N - m + 1)
            pattern = bin2dec(num2str(data(i:i + m - 1,j)')) + 1; % 获取当前 m 位模式
            observed_freq(pattern,j) = observed_freq(pattern,j) + 1; % 更新观察到的频率
        end
    end
    % 计算期望频率
    expected_freq = N / num_patterns; % 理想情况下的期望频率

    % 计算卡方统计量
    chi_squared = sum(((observed_freq - expected_freq).^2) / expected_freq);

    % 计算自由度
    df = num_patterns - 1; % 自由度为 (2^m - 1)

    % 计算 p 值
    p_value = 1 - chi2cdf(chi_squared, df); % 卡方分布的累积分布函数
    p_value_dis(2,9) = mean(p_value);

    save('figure_data\p_value_dis.mat','p_value_dis');
else
    load('figure_data\p_value_dis.mat');
end
%% NIST 显示
figure;
x = 1:9;
scatter(x, p_value_dis(1,:), 'filled');
hold on;
scatter(x, p_value_dis(2,:), 'filled');
xlabel('test number');
ylabel('p value');
plot([0,10],[0.01,0.01])
ylim([0.001, 1]); % 设置 y 轴范围
xlim([0.5, 9.5]); % 设置 y 轴范围
set(gca, 'YScale', 'log')

%% 编码自由度
% 输出编码自由度 (熵)
disp(['编码自由度 (DOF): ', num2str(mean(image_entropy_different))]);
disp(['信息编码容量: ', num2str(floor(radial_res*angular_res*2*mean(image_entropy_different)))]);






%% train
if first_time
    image_files_train = dir(fullfile('dataset\train\', '**','*.jpg'));
    train_data = false(radial_res*angular_res*2,length(image_files_train));
    train_labels = zeros(6,length(image_files_train));
    for ii = 1:length(image_files_train)
            % 获取图片文件的完整路径
            filenamestr = fullfile(image_files_train(ii).folder, image_files_train(ii).name);
            eyeimage_filename=filenamestr;
            %filenamestr='dataset\000\S6000S00.jpg';
            im = imread(filenamestr);
            I2 = imread(filenamestr);

%             eI=edge(I2,'canny', 0.2);

            % 利用hough变换找到图像中的一个圆
            %[y0detect,x0detect,Accumulator] = houghcircle(eI,45,4);

            dingwei();
            % with these settings a 9600 bit iris template is created
            guiyihua();

            tezhengtiqu();
            train_labels(str2num(filenamestr(end-8:end-7)),ii) = 1;
            train_data(:,ii) = template(:);%.*mask(:);
    end

    image_files_test = dir(fullfile('dataset\classify\', '**','*.jpg'));
    test_data = false(radial_res*angular_res*2,length(image_files_test));
    for ii = 1:length(image_files_test)
            % 获取图片文件的完整路径
            filenamestr = fullfile(image_files_test(ii).folder, image_files_test(ii).name);
            eyeimage_filename=filenamestr;
            %filenamestr='dataset\000\S6000S00.jpg';
            im = imread(filenamestr);
            I2 = imread(filenamestr);

%             eI=edge(I2,'canny', 0.2);

            % 利用hough变换找到图像中的一个圆
            [y0detect,x0detect,Accumulator] = houghcircle(eI,45,4);

            dingwei();
            % with these settings a 9600 bit iris template is created
            guiyihua();

            tezhengtiqu();
            test_data(:,ii) = template(:);%.*mask(:);
    end


    num_samples = size(train_data, 2);  % 样本数
    random_order = randperm(num_samples);  % 随机打乱样本索引
    train_data = train_data(:, random_order);
    train_labels = train_labels(:, random_order);

    % 定义BP神经网络
    hiddenLayerSize = 200; % 隐藏层神经元数量
    net = patternnet(hiddenLayerSize);

    % 禁用测试集，只保留训练集和验证集
    net.divideFcn = 'divideblock';  % 数据分割方式使用 'divideblock'
    net.divideParam.trainRatio = 0.8;  % 80% 数据用于训练
    net.divideParam.valRatio = 0.2;    % 20% 数据用于验证
    net.divideParam.testRatio = 0;     % 0% 数据用于测试（禁用测试集）

    % 设置训练参数
    net.trainParam.epochs = 100;  % 训练轮数
    net.trainParam.lr = 0.01;     % 学习率
    net.trainParam.min_grad = 0;     % 学习率
    % 禁用提前停止
    net.trainParam.max_fail = 100;  % 默认验证失败后停止训练，增加最大失败次数，避免过早停止


    % 训练神经网络
    [net, tr] = train(net, train_data, train_labels);

    % 训练完成后保存网络
    save('figure_data\trained_mnist_bp.mat', 'net', 'tr','test_data');
else
    load('figure_data\trained_mnist_bp.mat');
end
trainAcc = 1 - tr.perf;   % 训练集的准确率
valAcc = 1 - tr.vperf;    % 验证集的准确率

% 绘制训练和验证集的准确率曲线
figure;
epochs = 1:tr.num_epochs;  % 获取训练过程中有效的 epoch 数

plot(epochs, trainAcc(1:tr.num_epochs), '-o', 'LineWidth', 2);
hold on;
plot(epochs, valAcc(1:tr.num_epochs), '-x', 'LineWidth', 2);
xlabel('Epochs');
ylabel('Accuracy');
title('Training and Validation Accuracy');
legend('Training Accuracy', 'Validation Accuracy');
grid on;
hold off;

test_results = net( test_data);
test_results(test_results<0.9999) = test_results(test_results<0.9999)*0.1;

% 绘制三维柱状图
figure;
b = bar3(test_results(:,7:12),0.6);
colormap(jet);  % 可以尝试其他颜色映射，例如 hot, parula, etc.
caxis([0 1]);
% 为每个条形设置颜色
for k = 1:length(b)
    zdata = b(k).ZData;  % 获取 Z 轴数据，即条形图的高度
    b(k).CData = zdata;  % 将 Z 轴数据映射为颜色数据
    b(k).FaceColor = 'interp';  % 使用插值颜色，确保颜色平滑过渡
end
xticks([1 2 3 4 5 6]);  % 设置 X 轴刻度的位置
yticks([1 2 3 4 5 6]);  % 设置 Y 轴刻度的位置
yticklabels({'A1', 'A2', 'A3', 'A4', 'A5', 'A6'});  % X 轴刻度对应的标签
xticklabels({'B1', 'B2', 'B3', 'B4', 'B5', 'B6'});  % X 轴刻度对应的标签
zlim([0,1]);
xlim([0.5,6.5]);
ylim([0.5,6.5]);
title('Identical rate of genuine labels');
figure;
b=bar3(test_results(:,13:18),0.6);
colormap(jet);  % 可以尝试其他颜色映射，例如 hot, parula, etc.
caxis([0 1]);
% 为每个条形设置颜色
for k = 1:length(b)
    zdata = b(k).ZData;  % 获取 Z 轴数据，即条形图的高度
    b(k).CData = zdata;  % 将 Z 轴数据映射为颜色数据
    b(k).FaceColor = 'interp';  % 使用插值颜色，确保颜色平滑过渡   
end
xticks([1 2 3 4 5 6]);  % 设置 X 轴刻度的位置
yticks([1 2 3 4 5 6]);  % 设置 Y 轴刻度的位置
yticklabels({'A1', 'A2', 'A3', 'A4', 'A5', 'A6'});  % X 轴刻度对应的标签
xticklabels({'C1', 'C2', 'C3', 'C4', 'C5', 'C6'});  % X 轴刻度对应的标签
zlim([0,1]);
xlim([0.5,6.5]);
ylim([0.5,6.5]);
title('Identical rate of fake labels');