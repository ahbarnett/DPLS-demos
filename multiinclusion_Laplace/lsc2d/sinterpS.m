function A = sinterpS(n, y, L)
% SINTERPS interpolates the function f given at regular points, at any set
% of points y. (NUFFT could be used to accelerate this computation.)

c = 2*pi/L;

A = zeros(length(y), n);
for j = 1:n
    f = zeros(n,1); f(j) =1;
    fhat = fftshift(fft(f)/n);    
    for k=1:n
        A(:, j) = A(:, j) + fhat(k)*exp(1i*c*(-n/2+k-1)*y);
    end
end
    
A = real(A);

