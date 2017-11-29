function xout=Itp(x,p)
%
% This function interpolates a function of the form f(t)=q(t) at the
% points p on the interval [0,2*pi]. q(t) must be a periodic function.
%
% Example:
%   t=linspace(0,2*pi,33)';
%   s=2*pi*rand(100,1);
%   Itp(sin(t),s)-sin(s)
%

m=size(x,1)-1;
X=fft(x(1:m,:));
v=exp(1i*p);
z=zeros(size(p,1),size(x,2));
v0=v.^(-floor(m/2));
for k=1:m
    z=z+v0*X(mod(ceil(m/2)+k-1,m)+1,:)/m;
    v0=v0.*v;
end
xout=real(z);