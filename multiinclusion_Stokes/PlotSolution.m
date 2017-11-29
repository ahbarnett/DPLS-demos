function PlotSolution(FileName)

load(FileName)
s.x=zeros(sum(N),1);
s.len=N;
ls=0;
for k=1:length(N)
    s.x(ls+1:ls+s.len(k))=I{k}.x;
    ls=ls+s.len(k);
end
s=Quad(s);

uc.nei = 1; uc.d = 2*pi;
[x,w] = gauss(n); x = (1+x)/2; w = w'/2; % quadr on [0,1]
U.x=uc.d*x-uc.d/2+1i*uc.d/2; U.w=uc.d*w; U.nx=0*U.x-1i;
D=U; D.x=U.x-1i*uc.d;
H = uc.d; L.x = -uc.d/2+ 1i*H*x-uc.d/2*1i; L.nx = 0*L.x+1; L.w = H*w; % left side
R = L; R.x = L.x+uc.d; % right side


n=128;
[X,Y]=meshgrid(-pi:2*pi/(n-1):pi);

t.x=reshape(X,n^2,1)+1i*reshape(Y,n^2,1);

node=[];
for i=-1:1
    for j=-1:1
        node=[node;[real(s.x)+2*pi*i,imag(s.x)+2*pi*j]];
    end
end
edge=0*node;
ls=0;
for j=1:9
    for k=1:length(s.len)
        edge(ls+1:ls+s.len(k),:)=[(ls+1:ls+s.len(k))',[ls+2:ls+s.len(k),ls+1]'];
        ls=ls+s.len(k);
    end
end
in=inpoly([real(t.x),imag(t.x)],node,edge);
tt.x=t.x(~in);

if ~exist('uu')||size(uu,1)~=2*size(tt.x,1)
    uu=EvalSolution(FileName,tt);
    save(FileName,'uu','-append')
end

u=zeros(2*n^2,1);
u([~in;~in])=uu;
ux=u(1:end/2);
uy=u(end/2+1:end);
ux=reshape(ux,n,n);
uy=reshape(uy,n,n);
uu=sqrt(ux.^2+uy.^2);
figure
pcolor(X,Y,ux)
hold on
colormap(jet(1000))
axis equal
axis([-pi pi -pi pi])
shading flat
minu=min(min(uu));
maxu=max(max(uu));

ls=0;
for k=1:s.M
    for i=-1:1
        for j=-1:1
            fill(real(s.x([ls+1:ls+s.len(k),ls+1]))+2*pi*i,imag(s.x([ls+1:ls+s.len(k),ls+1]))+2*pi*j,'w','LineWidth',0.1)
        end
    end
   ls=ls+s.len(k);
end
caxis([minu,maxu])
title('u_x')

figure
pcolor(X,Y,uy)
hold on
colormap(jet(1000))
shading flat
axis equal
axis([-pi pi -pi pi])
ls=0;
for k=1:s.M
    for i=-1:1
        for j=-1:1
            fill(real(s.x([ls+1:ls+s.len(k),ls+1]))+2*pi*i,imag(s.x([ls+1:ls+s.len(k),ls+1]))+2*pi*j,'w','LineWidth',0.1)
        end
    end
    ls=ls+s.len(k);
end
caxis([minu,maxu])
title('u_y')

figure
pcolor(X,Y,uu)
hold on
colormap(jet(1000))
shading flat
axis equal
axis([-pi pi -pi pi])
ls=0;
for k=1:s.M
    for i=-1:1
        for j=-1:1
            fill(real(s.x([ls+1:ls+s.len(k),ls+1]))+2*pi*i,imag(s.x([ls+1:ls+s.len(k),ls+1]))+2*pi*j,'w','LineWidth',0.1)
        end
    end
    ls=ls+s.len(k);
end

caxis([minu,maxu])
title('|u|_2')
axis off

figure
ls=0;
for k=1:s.M
    for i=-1:1
        for j=-1:1
            fill(real(s.x([ls+1:ls+s.len(k),ls+1]))+2*pi*i,imag(s.x([ls+1:ls+s.len(k),ls+1]))+2*pi*j,'w','LineWidth',0.1)
            hold on
        end
    end
    ls=ls+s.len(k);
end
colormap(jet(1000))
axis equal
axis([-pi pi -pi pi])
axis off