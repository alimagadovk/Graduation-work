function h = Baterwort_filt(M,N,D0,n,x,y)
h = ones(M,N);
h1 = 1 - Baterwort_peak(M,N,D0,n,x,y);
h2 = 1 - Baterwort_peak(M,N,D0,n,-x,-y);
h = 1 - h .* h1 .* h2;
end