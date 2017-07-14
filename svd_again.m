function [ I2 ] = svd_again( I2,patch_size,L, step_patch,w,L_select_width,L_select_height,m,I,t,I_origin)
% patch_size:提取的patch的大小
% L：patch_group的大小
% step_patch：取patch时的step
% w：存放每个像素的权重矩阵
% L_select_width、L_select_height：在L_select_width * L_select_height范围内选取相似的L个patch
% m:每个patch的长与宽
% I：原始有噪声图片
% t：噪声的估计方差

[height, width] = size(I2);
[I2_sort,index] = sort(reshape(I2',1,width * height));                      % 对去噪图像的像素值进行升序排列
ind = 1;
while(1)
    if I2_sort(ind) < 0
        ind = ind + 1;
    else
        break;
    end
end
while(1)
    if I2_sort(ind) ~= 0
        break;
    else                                                                    % 存在像素值为0的点，则需要对其进行单独的svd处理
        x = floor(index(ind) / width) + 1;                                       % 计算得到该像素值的坐标位置
        y = index(ind) - width * (x - 1);
        if y == 0
            x = x - 1;
            y = width;
        end
        %% 得到范围框左上角点的坐标
        temp_L_height_index = x - round(L_select_height / 3);               % 计算范围框左上角点的纵坐标和横坐标（尽量使点x，y在中间）
        temp_L_width_index = y - round(L_select_width / 3);
        if temp_L_height_index < 1                                          % 若左上角点的纵坐标越界
            temp_L_height_index = 1;
        elseif temp_L_height_index + L_select_height - 1 > height           % 若左上角点的纵坐标造成范围框越界
            temp_L_height_index = height - L_select_height + 1;
        end
        if temp_L_width_index < 1                                          % 若左上角点的横坐标越界
            temp_L_width_index = 1;
        elseif temp_L_width_index + L_select_width - 1 > width           % 若左上角点的横坐标造成范围框越界
            temp_L_width_index = width - L_select_width + 1;
        end
        %% 得到包含未赋值像素点的中心patch的值
        temp_p_height_index = x - round(m / 3);               % 计算范围框左上角点的纵坐标和横坐标（尽量使点x，y在中间）
        temp_p_width_index = y - round(m / 3);
        if temp_p_height_index < 1                                          % 若左上角点的纵坐标越界
            temp_p_height_index = 1;
        elseif temp_p_height_index + m - 1 > height           % 若左上角点的纵坐标造成范围框越界
            temp_p_height_index = height - m + 1;
        end
        if temp_p_width_index < 1                                          % 若左上角点的横坐标越界
            temp_p_width_index = 1;
        elseif temp_p_width_index + m - 1 > width           % 若左上角点的横坐标造成范围框越界
            temp_p_width_index = width - m + 1;
        end
        center = reshape(I_origin(temp_p_height_index:temp_p_height_index + m - 1,temp_p_width_index:temp_p_width_index + m - 1)',1,m * m);             %%%%%
        [I2,w] = step1_2_3( temp_L_height_index,temp_L_width_index, patch_size,L,I2,step_patch,w,L_select_width,L_select_height,m,I,t,2,I_origin,center);
        %% 递归调用，直到不存在零像素点
        I2 = svd_again( I2,patch_size,L, step_patch,w,L_select_width,L_select_height,m,I,t,I_origin);
        break;
    end
end

end

