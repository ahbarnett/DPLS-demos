function L=ArcLength(x,y,l)
% This function computes the arc length of a vesicle with parameterization
% X(t) and Y(t). t must be define on an interval [0,l), where l is the
% width of the interval.
%
% Example:
%   t=linspace(0,2*pi,65)';
%   x=cos(t);
%   y=sin(t);
%   VWM_ArcLength(x,y,2*pi)
%

m=size(x,1);
sa=sqrt(Dtp(x).^2 + Dtp(y).^2);
L=l/(m-1)*sum(sa(1:m-1,:))';