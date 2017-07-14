function [ I2 ] = func_svd_lra4( I , t ,I_origin)

if t < 20
    patch_size = 9 * 9;                                                     % 根据噪声，确定提取patch时m的大小
elseif t < 40 && t >= 20
    patch_size = 10 * 10;
elseif t >= 40
    patch_size = 11 * 11;
end
m = sqrt(patch_size);
[height, width] = size(I);
I2 = zeros(size(I));                                                        % 存放新得到的图像
w = zeros(size(I));                                                         % 存放每个像素的权重矩阵
L = 85;
% L = 25;                                                                     % patch_group的大小
L_select_width = round(sqrt(L) * 4.1);                                      % 在L_select_width*L_select_height范围内选取相似的L个patch
L_select_height = round(sqrt(L) * 4.1);
% L_select_width = round(width / 13.5);                                      % 在L_select_width*L_select_height范围内选取相似的L个patch
% L_select_height = round(height / 13.5);
step_patch = 1;                                                             % 取patch时的step
step_L_select_width = round(floor(m) * 0.9);                                                    % 取L时范围框的横向移动步长
step_L_select_height = round(floor(m) * 0.9);                                                   % 取L时范围框的纵向移动步长
L_width_index = [1:step_L_select_width:width - L_select_width, width - L_select_width + 1];             % 计算范围框左上角横向的坐标
L_height_index = [1:step_L_select_height:height - L_select_height, height - L_select_height + 1];       % 计算范围框左上角纵向的坐标
for i = 1:length(L_height_index)
    for j = 1:length(L_width_index)
        [I2,w] = step1_2_3( L_height_index(i),L_width_index(j), patch_size,L,I2,step_patch,w,L_select_width,L_select_height,m,I,t,1 ,I_origin);
%         fprintf('%d,%d\n',i,j);
    end
end
imwrite(I2 / 255,'第一次去噪部分.jpg');
I2 = svd_again( I2,patch_size,L, step_patch,w,L_select_width,L_select_height,m,I,t,I_origin);
end