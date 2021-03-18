%% Тестирование методов фильтрации нормального белого шума
% накладываем нормальный белый шум
clc
clear
close all

im = imread('Lena.png');
figure
imshow(im)
title('Исходное изображение')
im = double(im);
pos_noise = imread('positive_noise.png');
neg_noise = imread('negative_noise.png');
noise = double(pos_noise) - double(neg_noise);
im_n = im + noise;
figure
imshow(uint8(im_n))
title('Зашумлённое изображение')

SNR_1 = snr(im,im_n - im);
MSE_1 = sum(sum((im - im_n).^2)) / (size(im,1) * size(im,2));
PSNR_1 = 20*log(255/sqrt(MSE_1))/log(10);
SSIM_1 = ssim(uint8(im),uint8(im_n));
%% тестирование частотных методов фильтрации нормального белого шума
clc

f = @(x)((cos(pi .* x) + 1) / 2).^1;
mask_size = 40;
ones_size = 8;
vect_down_len = (mask_size - ones_size) / 2;
mask = SquareMask(mask_size,ones_size,f,0,0);


im_frame = padarray(im_n, [vect_down_len vect_down_len], 'symmetric', 'both');
n_frame = padarray(noise, [vect_down_len vect_down_len], 'symmetric', 'both');

% Выбираем фильтр
im_2 = Filtr_Winn(im_frame,noise,mask,vect_down_len,2);
%im_2 = Hard_Filt(im_frame,mask,vect_down_len,2,0.16);
%im_2 = Soft_Filt(im_frame,mask,vect_down_len,2,0.076);


i1 = vect_down_len + 1;
i2 = size(im_frame,1) - vect_down_len;
j1 = vect_down_len + 1;
j2 = size(im_frame,2) - vect_down_len;
im_2 = im_2(i1:i2,j1:j2);

figure
imshow(uint8(im_2))
title('Отфильтрованное изображение')



SNR_2 = snr(im,im_2 - im);
MSE_2 = sum(sum((im - im_2).^2)) / (size(im,1) * size(im,2));
PSNR_2 = 20*log(255/sqrt(MSE_2))/log(10);
SSIM_2 = ssim(uint8(im),uint8(im_2));
disp('SNR (' + string(SNR_1) + '/' + string(SNR_2) + ')')
disp('PSNR (' + string(PSNR_1) + '/' + string(PSNR_2) + ')')
disp('SSIM (' + string(SSIM_1) + '/' + string(SSIM_2) + ')')
%% тестирование пространственных методов фильтрации нормального белого шума
clc
N = 3;
% Выбираем фильтр
im_res = spfilt(im_n, 'amean', N, N);
%im_res = spfilt(im_n, 'midpoint', N, N);

% D0 = 0.9;
% f = fspecial('gaussian', [N N], D0);
% im_res = imfilter(im_n,f,'replicate');

%im_res = bilateral_filter1(im_n,3,1.2,80);

figure
imshow(uint8(im_res))
title('Отфильтрованное изображение')

SNR_2 = snr(im,im_res - im);
MSE_2 = sum(sum((im - im_res).^2)) / (size(im,1) * size(im,2));
PSNR_2 = 20*log(255/sqrt(MSE_2))/log(10);
SSIM_2 = ssim(uint8(im),uint8(im_res));
disp('SNR (' + string(SNR_1) + '/' + string(SNR_2) + ')')
disp('PSNR (' + string(PSNR_1) + '/' + string(PSNR_2) + ')')
disp('SSIM (' + string(SSIM_1) + '/' + string(SSIM_2) + ')')
%% Тестирование методов фильтрации периодической помехи
% накладываем периодическую помеху на изображение
clc
clear
close all
im = imread('Lena.png');
im = double(im);
F1 = 100.4;
F2 = 100.2;
C = [F1 F2]; % координаты частот, составляющих гарм. помехи на спектре относительно центра изображения
A(1:size(C,1)) = [90]; % амплтитуды составл. гарм. помехи
[r, R, S] = MyGarmNoise(512, 512, C, A);
im_n = im + r;

figure
imshow(uint8(im))
title('Исходное изображение')

figure
imshow(uint8(im_n))
title('Изображение с наложенной гарм. помехой')

SNR_1 = snr(im,im_n - im);
MSE_1 = sum(sum((im - im_n).^2)) / (size(im,1) * size(im,2));
PSNR_1 = 20*log(255/sqrt(MSE_1))/log(10);
SSIM_1 = ssim(uint8(im),uint8(im_n));
%% тестирование реализации фильтра Винера для подавления периодической помехи
clc
im_2 = Filtr_Winn_garm(im_n,R);

figure
imshow(uint8(im_2))
title('Отфильтрованное изображение')


SNR_2 = snr(im,im_2 - im);
MSE_2 = sum(sum((im - im_2).^2)) / (size(im,1) * size(im,2));
PSNR_2 = 20*log(255/sqrt(MSE_2))/log(10);
SSIM_2 = ssim(uint8(im),uint8(im_2));
disp('SNR (' + string(SNR_1) + '/' + string(SNR_2) + ')')
disp('PSNR (' + string(PSNR_1) + '/' + string(PSNR_2) + ')')
disp('SSIM (' + string(SSIM_1) + '/' + string(SSIM_2) + ')')
%% тестирование метода оптимальной узкополосной фильтрации
clc
M = 512;
N = 512;

% Выбираем функцию для выделения шумовой составляющей
D0 = 25;
n = 11;
h = Baterwort_filt(M,N,D0,n,F1,F2);

% D0 = 22;
% h = Gauss_filt(M,N,D0,F1,F2);


spectr = fftshift(fft2(im_n));

spectr_n = spectr .* h;


noise = ifft2(ifftshift(spectr_n));

figure
imshow(uint8(noise))
title('Полученное изображение помехи')

PSF = fspecial('gaussian',10,21);

im_n = real(im_n);
noise = real(noise);

im_n = edgetaper(im_n,PSF);
noise = edgetaper(noise,PSF);

im_2 = Opt_filtr(im_n,noise,8);

figure
imshow(uint8(im_2))


SNR_2 = snr(im,im_2 - im);
MSE_2 = sum(sum((im - im_2).^2)) / (size(im,1) * size(im,2));
PSNR_2 = 20*log(255/sqrt(MSE_2))/log(10);
SSIM_2 = ssim(uint8(im),uint8(im_2));
disp('SNR (' + string(SNR_1) + '/' + string(SNR_2) + ')')
disp('PSNR (' + string(PSNR_1) + '/' + string(PSNR_2) + ')')
disp('SSIM (' + string(SSIM_1) + '/' + string(SSIM_2) + ')')