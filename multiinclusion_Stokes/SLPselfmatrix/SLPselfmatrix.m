function S=SLPselfmatrix(s)
%
% This function generates the single-layer self interaction matrix using
% the Hybrid Gauss-Trapezoidal rule.
%
% Example:
%   t=linspace(0,2*pi,65)'; % Note that both endpoints are included.
%   s=[cos(t)+sin(2*t)/4,sin(t)/2];
%   S=SLPselfmatrix(s);

N=size(s,1); t=linspace(0,2*pi,N)'; C=tril(ones(N-1)); C=logical([C;~C]);
Y1=s(2:end,1)*ones(1,N-1); Y1=[Y1;Y1]; Y1=reshape(Y1(C),[N-1,N-1]).';
Y2=s(2:end,2)*ones(1,N-1); Y2=[Y2;Y2]; Y2=reshape(Y2(C),[N-1,N-1]).';
Y1=[Y1,Y1(:,1)]; Y1=[Y1(end,:);Y1]; Y2=[Y2,Y2(:,1)]; Y2=[Y2(end,:);Y2];
Y=[Y1,Y2]; Y=[Y,Dtp(Y)]; G=@(x,y)SLPG(y,s);

% Alpert quadrature is done here.
S=Alpert(G,t,Y,1);

S11=S(1:end/3,:); S11(:,end)=S11(:,1)+S11(:,end); S11=S11(2:end,2:end);
S12=S(end/3+1:2*end/3,:); S12(:,end)=S12(:,1)+S12(:,end);
S12=S12(2:end,2:end); S22=S(2*end/3+1:end,:);
S22(:,end)=S22(:,1)+S22(:,end); S22=S22(2:end,2:end); C=C(:,end:-1:1);
S11=[S11,S11].'; S11=reshape(S11(C),[N-1,N-1]).';
S12=[S12,S12].'; S12=reshape(S12(C),[N-1,N-1]).';
S22=[S22,S22].'; S22=reshape(S22(C),[N-1,N-1]).';
S=[S11,S12;S12,S22];
end

function z=SLPG(y,s)
%
% The single layer kernel is defined here.
%

N=size(y,1); Y1=y(:,1:end/4); Y2=y(:,end/4+1:end/2);
DY1=y(:,end/2+1:3*end/4); DY2=y(:,3*end/4+1:end);
R1=ones(N,1)*s(:,1)'-Y1; R2=ones(N,1)*s(:,2)'-Y2;
RHO2=R1.^2+R2.^2; LOGRHO=-log(sqrt(RHO2)); SA=sqrt(DY1.^2+DY2.^2);
z11=1/4/pi*(LOGRHO+(R1.*R1)./RHO2).*SA;
z12=1/4/pi*((R1.*R2)./RHO2).*SA;
z22=1/4/pi*(LOGRHO+(R2.*R2)./RHO2).*SA;
z=[z11,z12,z22];
end