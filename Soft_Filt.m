function [im_res] = Soft_Filt(im_frame,mask,vect_down_len,k_m,lamb)
i1 = 1 + vect_down_len;
i2 = size(im_frame,1) - vect_down_len;
j1 = 1 + vect_down_len;
j2 = size(im_frame,2) - vect_down_len;
im_n = im_frame(i1:i2,j1:j2);
im_res = zeros(size(im_frame));
step = size(mask,1) - vect_down_len;
lamb = lamb*(k_m*size(mask,1))^2;
for i = 0:step:size(im_n,1)
    for j = 0:step:size(im_n,2)
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
        for k = 1:size(sp_window,1)
            mas1 = find(abs(sp_window(k,:)) < lamb);
            mas2 = find(abs(sp_window(k,:)) > lamb);
            sp_window(k,mas1(:)) = 0;
            for l = 1:length(mas2)
                z = sp_window(k,mas2(l));
                koeff = (abs(z) - lamb)/(abs(z));
                sp_window(k,mas2(l)) = z.*koeff;
            end            
        end
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