function [ I2,w ] = step1_2_3( L_height_index,L_width_index, patch_size,L,I2,step_patch,w,L_select_width,L_select_height,m ,I,t,flag,I_origin,center )
% L_height_index：范围框左上角纵向的坐标
% L_width_index：范围框左上角横向的坐标
% patch_size:提取的patch的大小
% L：patch_group的大小
% step_patch：取patch时的step
% w：存放每个像素的权重矩阵
% L_select_width、L_select_height：在L_select_width * L_select_height范围内选取相似的L个patch
% m:每个patch的长与宽
% I：原始有噪声图片
% t：噪声的估计方差
global k_record_index_iter
global k_record
global k_record_index
%% step1——patch grouping: block matching method
L_index_begin = [L_height_index L_width_index];                                                     % 取L范围框左上角点的坐标
patch_height_index = [L_index_begin(1):step_patch:L_index_begin(1) + L_select_height - 1 - m, ...    % 范围框内patch左上角初始位置的索引
    L_index_begin(1) + L_select_height - m];
patch_width_index = [L_index_begin(2):step_patch:L_index_begin(2) + L_select_width - 1 - m, ...    % 范围框内patch左上角初始位置的索引
    L_index_begin(2) + L_select_width - m];
ind = 0;
number_patch = length(patch_width_index) * length(patch_height_index);
y = zeros(number_patch, patch_size);
y2 = zeros(number_patch, patch_size);
P = zeros(patch_size, L + 1);
for k1 = 1:m
    for k2 = 1:m
        ind = ind + 1;
        y(:,ind) = reshape(I(patch_height_index + k1 - 1,patch_width_index + k2 - 1)',1,number_patch);                % 得到由范围框内的像素构成的patch
        y2(:,ind) = reshape(I_origin(patch_height_index + k1 - 1,patch_width_index + k2 - 1)',1,number_patch);  
    end
end
% 使用欧几里得距离，对中心patch，在范围框内，寻找其L个最近邻
center_patch_num = 2;
for iii = 1:center_patch_num^2
    if flag == 1
        center_patch_x = floor((iii - 1) / center_patch_num) + 1;
        center_patch_y = iii - (center_patch_x - 1) * center_patch_num; 
        temp = round(length(patch_height_index) / (center_patch_num * 2) * (center_patch_x * 2 - 1) - 1) * length(patch_height_index) + ...
            round(length(patch_width_index) / (center_patch_num * 2) * (center_patch_y * 2 - 1));   % 取范围框中间的patch作为中心patch
        patch_center = y2(temp,:);                  %%%%%%
    else
        patch_center = center;
    end
    % 计算当前所选模块的噪声方差 
    % [~,B,C,D] = dwt2(I(L_index_begin(1):L_index_begin(1) + L_select_height - 1,...
    %     L_index_begin(2):L_index_begin(2) + L_select_width - 1),'db1');                                              % 二维小波变换
    % t_L = median(median(median(abs([B,C,D])))) / 0.6745 ;                   % 估计噪声的标准偏差
    % 
    % t_par = (t / t_L)^3;
    % L = round(L * t_par);
    % if L >= number_patch
    %     L = number_patch - 1;
    % end
    % fprintf('%d\n',L);
%     temp_dis_norm = (y2 - repmat((mean(y2'))',1,m * m))./ repmat((var(y2'))',1,m * m)...
%         - repmat((patch_center - repmat(mean(patch_center),1,m * m))./repmat(var(patch_center),1,m*m),number_patch,1);      % 先归一化之后再计算欧几里得距离，与原文不同
    y2 = y2(:,:) - repmat(mean(y2,2),1,patch_size) * 0.2;
    patch_center = patch_center - repmat(mean(patch_center),1,patch_size) * 0.2;    
    temp_dis = y2(:,:) -  repmat(patch_center,number_patch,1);               % 计算其他的patch到中心patch之间的欧几里得距离                  %%%%
    temp_dis2 = sum(temp_dis.^2,2);
    [~,index] = sort(temp_dis2(:));                                         % 按照欧几里得距离进行排序
    P(:,1:L + 1) = y(index(1:L + 1),:)';                                    % 按距离从近到远取出L个patch
    %% step2——SVD_based denoising
    P_ave = mean(P');
    P_SVD = zeros(size(P));
    P = P - repmat(P_ave',1,L + 1);
    [U,S,V] = svd(P);                                                       % 进行SVD变换
    [M,N] = size(P);
    [n1,n2] = size(S);
    sum_k_1 = 0;
    k = 0;
    for j2 = min(n1,n2) : -1 : 2
        sum_k_1 = sum_k_1 + S(j2,j2)^2;                                     % 对于通过svd变换得到的特征值，从后往前求取平方和
        sum_k = sum_k_1 + S(j2 - 1,j2 - 1)^2;
        if sum_k_1 <= t^2 * (L + 1) * patch_size && sum_k >= t^2 * (L + 1) * patch_size             % 当满足论文所给条件之后
            k = j2 - 1;                                                     % 获得应取的维数
            break;
        end
    end
    S2 = zeros(k,k);
    U2 = zeros(M,k);
    V2 = zeros(N,k);
    S2(:,:) = S(1:k,1:k);                                                   % 按照维数r重构出矩阵
    U2(:,:) = U(:,1:k);
    V2(:,:) = V(:,1:k);
    P_SVD(:,:) = U2 * S2 * V2' + repmat(P_ave',1,L + 1);                                                 % 将变换后的patch group保存到矩阵P_SVD 
%     fprintf('%d,%d\n',k,L);
    P_SVD_error = P_SVD - (P + repmat(P_ave',1,L + 1));
    P_SVD_error_var = var(P_SVD_error);
    P_SVD_select = (P_SVD_error_var <= (t^2 * 0.8)); 
    k_record_index = k_record_index + 1;
    k_record(k_record_index,k_record_index_iter) = k;
%     fprintf('%d\n',k);
    %% step3——aggregation
    % aggregate pixel
    w_temp = 0;
    if k < L + 1
        w_temp = 1 - k/(L + 1);                                             % 计算该像素的权重
    elseif k == L + 1
        w_temp = 1/(L + 1);
    end
    for ii = 1 : L + 1
        %if P_SVD_select(ii) == 1 || flag == 2
            index_height = floor(index(ii) ./ length(patch_width_index)) + 1;
            index_width = index(ii) - length(patch_width_index) * (index_height - 1);
            if index_width == 0
                index_width = length(patch_width_index);
                index_height = index_height - 1;
            end
            X = patch_height_index(index_height);
            Y = patch_width_index(index_width);
            temp = w(X:X + m - 1,Y:Y + m - 1) + w_temp .* ones(m,m);
            I2(X:X + m - 1, Y: Y + m - 1) = ((I2(X:X + m - 1, Y: Y + m - 1) .* w(X:X + m - 1,Y:Y + m - 1) ...
                    + reshape(P_SVD(:,ii)',m,m)' * w_temp)) ./ temp;
            w(X:X + m - 1,Y:Y + m - 1) = temp;
       % end
    end
    if flag ~= 1
        break;
    end
end
end

