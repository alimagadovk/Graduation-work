function mask = SquareMask(size_mask,size_ones,func,vy,flag)
if (size_ones > size_mask)
    error('Длина единичной части должна быть меньше длины маски!')
end
if (mod(size_ones, 2) ~= mod(size_mask, 2))
    error('Чётность длин маски и единичной части должна совпадать!')
end
if (flag)
    vect_down = vy;
    i1 = (size_mask - size_ones)/2;
    i2 = (size_mask + size_ones)/2;
    if (i1 ~= length(vect_down))
        error('Неверно задана длина вектора спуска!')
    end
else
    vect_down_len = (size_mask - size_ones)/2;
    i1 = vect_down_len;
    i2 = size_mask - vect_down_len;
    step = 1 / (vect_down_len - 1);
    vx = 0:step:1;
    vect_down = func(vx);
end
lent_mask = zeros(1,size_mask);
lent_mask(i1 + 1:i2) = 1;
vect_down = sort(vect_down);
lent_mask(1:i1) = vect_down(:);
lent_mask(i2 + 1:end) = vect_down(end:-1:1);
mask = lent_mask' * lent_mask;
f_test = lent_mask(1: i1) + lent_mask(i2 + 1: end);
f_test = sum(f_test)/length(f_test);
if (f_test ~= 1)
    error('Сумма перекрывающихся частей не равна 1!')
end
end