function PlotSolution(FileName)
%
% This function plots the solution stored in FileName.

load(FileName)
if ~exist('P')
   P=M;
   M=n;
   clear n
end
nei = 1;  

s.x=zeros(sum(N),1);
s.len=N;
ls=0;
for k=1:length(N)
    s.x(ls+1:ls+s.len(k))=I{k}.x/2/pi;
    ls=ls+s.len(k);
end
s=Quad(s);

U.e1 = 1; U.e2 = 1i; U.nei = 0; 
[x,w] = gauss(M); w=w/2;  % walls
L.x = (-U.e1 + U.e2*x)/2; L.nx = (-1i*U.e2)/abs(U.e2) + 0*L.x; L.w=w*abs(U.e2);
R = L; R.x = L.x + U.e1;
B.x = (-U.e2 + U.e1*x)/2; B.nx = (1i*U.e1)/abs(U.e1) + 0*B.x; B.w=w*abs(U.e1);
T = B; T.x = B.x + U.e2;
Rp = 1.4; %P=128;
p.x = Rp * exp(1i*(1:P)'/P*2*pi); p = quadr(p); % proxy pts
U.L=L;
U.R=R;
U.T=T;
U.B=B;



n=256;
[X,Y]=meshgrid(-0.5:1/(n-1):0.5);

t.x=reshape(X,n^2,1)+1i*reshape(Y,n^2,1);

%if ~exist('u','var')||size(u,1)~=size(t.x,1)   % ahb considers dangerous.
    u=zeros(n^2,1);
    ux=u;
    uy=u;
    ls=0;
    while ls<n^2
        lsmax=min(ls+2e5,n^2);
        tt.x=t.x(ls+1:lsmax);
        [uu,uux,uuy]=At_operator(tt,s,Sig,U,nei);
        u(ls+1:lsmax)=u(ls+1:lsmax)+LapSLPmatrix(tt,p)*Psi+uu;
        tt.nx=0*tt.x+1;
        [~,UX]=LapSLPmatrix(tt,p);
        ux(ls+1:lsmax)=ux(ls+1:lsmax)+UX*Psi+uux;
        tt.nx=0*tt.x+1i;
        [~,UY]=LapSLPmatrix(tt,p);
        uy(ls+1:lsmax)=uy(ls+1:lsmax)+UY*Psi+uuy;
        ls=lsmax;
    end
%    save(FileName,'u','ux','uy','-append')
%end

u=reshape(u,n,n);
ux=reshape(ux,n,n);
uy=reshape(uy,n,n);
figure
pcolor(X,Y,u)
hold on
colormap(jet(1000))
axis equal
axis([-0.5 0.5 -0.5 0.5])
shading flat
minu=min(min(u));
maxu=max(max(u));
LCurv=linspace(minu,maxu,10*sqrt(length(s.len))/5);
LCurv=LCurv(2:end-1);
contour(X,Y,u,LCurv,'k','LineWidth',0.5)
colorbar; title('u')
set(gcf, 'Renderer', ' painters');
ls=0;
for k=1:s.M
    for i=-1:1
        for j=-1:1
            fill(real(s.x([ls+1:ls+s.len(k),ls+1]))+i,imag(s.x([ls+1:ls+s.len(k),ls+1]))+j,'w','LineWidth',0.1)
        end
    end
   ls=ls+s.len(k);
end

if 0
figure
ls=0;
for k=1:s.M
    for i=-1:1
        for j=-1:1
            fill(real(s.x([ls+1:ls+s.len(k),ls+1]))+i,imag(s.x([ls+1:ls+s.len(k),ls+1]))+j,'w','LineWidth',0.1)
            hold on
        end
    end
    ls=ls+s.len(k);
end
colormap(jet(1000))
axis equal
axis([-0.5 0.5 -0.5 0.5])
axis off
set(gca,'color','none')
end

figure
pcolor(X,Y,sqrt(ux.^2+uy.^2))
hold on
colormap(jet(1000))
shading flat
axis equal
axis([-0.5 0.5 -0.5 0.5])
ls=0;
for k=1:s.M
    for i=-1:1
        for j=-1:1
            fill(real(s.x([ls+1:ls+s.len(k),ls+1]))+i,imag(s.x([ls+1:ls+s.len(k),ls+1]))+j,'w','LineWidth',0.1)
        end
    end
    ls=ls+s.len(k);
end
uxmin=min(min(sqrt(ux.^2+uy.^2)));
uxmax=max(max(sqrt(ux.^2+uy.^2)));
%caxis([uxmin,(uxmax+uxmin)/8]);
colorbar; title('|grad u|')
end