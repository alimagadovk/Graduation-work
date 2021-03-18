function h = Baterwort_peak(M,N,D0,n,x,y)
if nargin == 4
    x = 0;
    y = 0;
end
D = @(u,v)sqrt((u - M/2 - 1 - x).^2 + (v - N/2 - 1 - y).^2);
H = @(u,v) 1./(1 + (D(u,v)./D0).^(2*n));
X = 1:M;
Y = 1:N;
[X,Y] = meshgrid(X,Y);
h = H(X,Y);
end