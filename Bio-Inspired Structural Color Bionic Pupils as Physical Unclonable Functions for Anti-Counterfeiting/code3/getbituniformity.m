function [entropy,bit_uniformity] = getbituniformity(template, mask)
    % 获取非噪声区域，即 mask 为 0 的区域
    valid_template = template(mask == 0);
    
    % 计算0和1的个数
    num_zeros = sum(valid_template == 0);
    num_ones = sum(valid_template == 1);
    
    % 非噪声区域的总元素个数
    total_elements = numel(valid_template);
    
    % 如果总元素个数为0，熵为0
    if total_elements == 0
        entropy = 0;
        return;
    end
    
    % 计算0和1的概率
    p0 = num_zeros / total_elements;
    p1 = num_ones / total_elements;
    
    % 避免log2(0)的情况
    if p0 == 0
        entropy0 = 0;
    else
        entropy0 = -p0 * log2(p0);
    end
    
    if p1 == 0
        entropy1 = 0;
    else
        entropy1 = -p1 * log2(p1);
    end
    
    % Shannon熵
    entropy = entropy0 + entropy1;
    bit_uniformity=p0;
end


