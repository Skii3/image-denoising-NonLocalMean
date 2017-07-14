% 论文：An efficient SVD_based method for image denoising
close all
clear
clc
clear global
global k_record
global k_record_index
global k_record_index_iter
k_record_index = 0;
k_record_index_iter = 1;
%% read the noise image
for noise_sigma = 10:10:10
    % I = imread('lena.jpg');                                           % 读取图片
    I = imread('fingerprint.png'); 
% I = imread('3.jpg'); 
    I_origin = double(I);
    temp = size(I);                                                    % 得到通道数，若为3通道则转化为1通道，若为1通道则不处理
    if length(temp) == 3 && temp(3) == 3
        I = rgb2gray(I);
    end
    imwrite(I,'原图像.jpg');                                            % 保存原始图像
    imshow(I),title('原始图像');                                        % 显示原始图像

    randn('seed', 0);                                                   % 生成随机种子
    noise = noise_sigma * randn(size(I));                                        % 通过随机种子得到噪声图像
    I = double(I);
    I = I + noise;                                                      % 得到加入噪声之后的图
    I = imresize(I, [512,512]);
     temp = size(I);   
%     I_origin = double(I);
    figure,imshow(I / 255),title('加高斯白噪声');                        % 保存加入噪声之后的图
    imwrite(I / 255,'加高斯白噪声.jpg');
    tic
    %% estimate noise standard devition by computing median absolute devition of the finest wavelet coefficients
    %% 进行高斯低通滤波
    f = fspecial('gaussian',[2 2],5);    
    I_origin_est = imfilter(I,f,'same');

    figure,imshow(I_origin_est,[]);

    [~,B,C,D] = dwt2(I,'db1');                                          % 二维小波变换
    t = median(median(median(abs([B,C,D])))) / 0.6745 ;                  % 估计噪声的标准偏差
%     scale = t_scale(t);
%     t = t * scale ;                                                     % 噪声估计的校准
    I_denoise = func_svd_lra4(I,t,I_origin_est);                        % 调用函数进行svd去噪，I_origin_est是patch group的基准对象
    tab1 = tabulate(k_record(:,1));
    figure,imshow(I_denoise,[]),title('第一次去噪图像');
    imwrite(I_denoise / 255,'第一次去噪.jpg');
    K = I_denoise;
    MSE = sum(sum((I_origin - K).^2)) / temp(1) / temp(2);
    PSNR1 = 20 * log10(255 / sqrt(MSE));
    %% back projection
    k_record_index_iter = 2;
    delta = 0.5;                                                          
    I_back = I_denoise + delta .* (I - I_denoise);                      % 获得新的有噪声的图像
    I_origin_est = imfilter(I_back,f,'same');                           % 进行滤波，简单去噪
    figure,imshow(I_origin_est,[]);
%     for gamma = 0.60:0.01:0.60
    gamma = 0.60;
        % t_2 = gamma * sqrt(abs(t^2 - sum(sum((I - I_back).^2))/(temp(1) * temp(2))));     
%         [~,B,C,D] = dwt2(I_back,'db1');                                     % 二维小波变换
%         t_2 = median(median(median(abs([B,C,D])))) / 0.6745 * 0.8;                % 估计新图像的噪声的标准偏差
    %     scale = t_scale2(t_2);
    %     t_2 = t_2 * scale;
        t_2 = gamma * sqrt(abs((noise_sigma^2 * temp(1) * temp(2) - sum(sum(I - I_back).^2)))/ temp(1) / temp(2));
        I_final = func_svd_lra4(I_back,t_2,I_origin_est);                   % 再一次调用函数进行svd去噪
        tab2 = tabulate(k_record(:,2));
        figure,imshow(I_final,[]),title('第二次去噪图像');
        imwrite(I_final / 255,'第二次去噪.jpg');                             
        figure,imshow(I - I_denoise,[]);                                    %显示原始图像与第一次去噪图像的差
        figure,imshow(I_final - I_back,[]);                                 %显示第二次去噪图像与新的有噪声图像之差
        toc
        %% 计算信噪比
        K = I_final;
        MSE = sum(sum((I_origin - K).^2)) / temp(1) / temp(2);
        PSNR = 20 * log10(255 / sqrt(MSE))
%         close all
%         gamma
%     end
end