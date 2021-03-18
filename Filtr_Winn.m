function im_res = Filtr_Winn(im_frame,noise,mask,vect_down_len,k_m)
N = size(im_frame,1) - 2*vect_down_len;
M = size(im_frame,2) - 2*vect_down_len;
sp_n = fftshift(fft2(noise));
Sn = size(mask,1)*size(mask,2)*var(sp_n(:))/(M*N);
im_res = zeros(size(im_frame));
step = size(mask,1) - vect_down_len;
for i = 0:step:N
    for j = 0:step:M
        i1 = i + 1;
        i2 = i + size(mask,1);
        j1 = j + 1;
        j2 = j + size(mask,2);
        im_window = double(im_frame(i1:i2,j1:j2)) .* mask;
        im_temp = zeros(k_m.*size(mask));
        i1 = (size(im_temp,1) - size(im_window,1))/2;
        i2 = (size(im_temp,1) + size(im_window,1))/2;
        j1 = (size(im_temp,2) - size(im_window,2))/2;
        j2 = (size(im_temp,2) + size(im_window,2))/2;
        im_temp(i1 + 1:i2,j1 + 1:j2) = im_window;
        sp_window = fftshift(fft2(im_temp));
        Sg = sp_window .* conj(sp_window);
        Sf = Sg - Sn;
        for k = 1:size(Sf,1)
            mas = find(Sf(k,:) < 0);
            for l = 1:length(mas)
                Sf(k,mas(l)) = 0;
            end
        end
        Filt_winn = Sf ./ Sg;
        mas = find(Sg == 0);
        Filt_winn(mas) = 0;
        Filt_winn(floor(size(sp_window,1)/2) + 1,floor(size(sp_window,2)/2) + 1) = 1;
        sp_window = Filt_winn .* sp_window;
        im_temp = ifft2(ifftshift(sp_window));
        im_window = im_temp(i1 + 1:i2,j1 + 1:j2);                
        i1 = i + 1;
        i2 = i + size(mask,1);
        j1 = j + 1;
        j2 = j + size(mask,2);
        im_res(i1:i2,j1:j2) = im_res(i1:i2,j1:j2) + im_window;
    end
end
end