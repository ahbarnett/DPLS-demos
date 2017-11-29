function z=Alpert(G,t,y,d)

m=length(t); if nargin==4&&d==1, L0=t(end)-t(1); L=L0/2; m0=m; m=(m+1)/2; 
else L=t(end)-t(1); end
% Alpert nodes and weights for 8th order integration.
v=[6.531815708567918e-3,9.086744584657729e-2,3.967966533375878e-1,...
    1.027856640525646e0,1.945288592909266e0,2.980147933889640e0,...
    3.998861349951123e0]';
u=[2.462194198995203e-2,1.701315866854178e-1,4.609256358650077e-1,...
    7.947291148621894e-1,1.008710414337933e0,1.036093649726216e0,...
    1.004787656533285e0]';
x=[2.087647422032129e-1,9.786087373714483e-1,1.989541386579751e0,3e0]';
w=[5.207988277246498e-1,9.535038018555888e-1,1.024871626402471e0,...
    1.000825744017291e0]';
g=11;n=round(m-9);h=L/(n+8);a=5;

if nargin==4&&d==1
    IT=Itp(eye(m0),[v*h;L-x*h]);
    tg=[v*h;L-x*h;t(a+1:a+n)];tg=[tg;L0-tg];
    if isempty(y), Geval=h*G(tg).'; else yud=y(end:-1:1,:);
    yg=[IT*y;y(a+1:a+n,:);IT*yud;yud(a+1:a+n,:)]; Geval=h*G(tg,yg).'; end
    W=(ones(size(Geval,1),1)*[u;w]');
    zL=(W.*Geval(:,1:g))*IT;
    zR=(W.*Geval(:,end/2+1:end/2+g))*IT;
    zL(:,a+1:a+n)=zL(:,a+1:a+n)+Geval(:,g+1:end/2);
    zR(:,a+1:a+n)=zR(:,a+1:a+n)+Geval(:,end/2+g+1:end);
    z=zL+zR(:,end:-1:1);
else
    IT=Itp(eye(m),[v*h;L-x*h]);
    tg=[v*h;L-x*h;t(a+1:a+n)]; if isempty(y), Geval=h*G(tg).'; else
    yg=[IT*y;y(a+1:a+n,:)]; Geval=h*G(tg,yg).'; end
    z=((ones(size(Geval,1),1)*[u;w]').*Geval(:,1:g))*IT;
    z(:,a+1:a+n)=z(:,a+1:a+n)+Geval(:,g+1:end);
end