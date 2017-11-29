function dxdt=Dtp(x,d)

persistent K
%
% This program computes the derivative of a function of the form 
% f(t)=a*t+p(t), where a is a constant and p(t) is a periodic function on 
% the interval [0,2*pi]. Both endpoints must be included.
%
% Example:
%   t=linspace(0,2*pi,33)';
%   y=sin(t);
%   dy=Dt(y);
%

[m,p]=size(x);
m=m-1;
if isempty(K)||size(K,1)~=m||size(K,2)~=p
    if (-1)^m==1
        K=1i*[(0:m/2-1)';0;(-m/2+1:-1)']*ones(1,p);
    else
        K=1i*[(1:(m-1)/2)';(-(m-1)/2:0)']*ones(1,p);
    end
end

if nargin==1
    dxdt=ifft(K.*fft(x(1:m,:)));
else
    dxdt=ifft(K.^d.*fft(x(1:m,:)));
end
dxdt=real([dxdt;dxdt(1,:)]);