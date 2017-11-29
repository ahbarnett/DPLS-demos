function [I,N]=GenerateIslands(M,n)

S=[];
R0=10*sqrt(10)/6.9/sqrt(M);
rate=1.05;
R=R0*rate;
while size(S,1)/n<M
    R=R/rate;
    MR=round(5/R^2);
    %S0=MakeSofR(R,MR,n);
    S0=MakeBarnettShape(R,MR,n);
    S0=IslandIntersect(S0,S,n);
    S=IslandIntersect(S,S0,n);
    S0=RecursiveIntersect(S0,n);
    S0=RecursiveIntersect(S0,n);
    S=[S;S0];
end
S=S(1:M*n,:);

ls=0;
I=cell(M,1);
for k=1:M
    I{k}.x=S(ls+1:ls+n,:)*[1;1i]-pi-1i*pi;
    area=Area(real(I{k}.x([1:end,1])),imag(I{k}.x([1:end,1])));
    xc=Area(real(I{k}.x([1:end,1])).^2/2,imag(I{k}.x([1:end,1])))/area;
    yc=Area(real(I{k}.x([1:end,1])).*imag(I{k}.x([1:end,1])),imag(I{k}.x([1:end,1])))/area;
    I{k}.x=0.97*(I{k}.x-xc-1i*yc)+xc+1i*yc;
    ls=ls+n;
end
N=n*ones(length(I),1);


ls=0;
for k=1:M
    for i=-1:1
        for j=-1:1
            plot(I{k}.x([1:end,1])+2*pi*i+1i*2*pi*j)
            hold on
        end
    end
    ls=ls+n;
end
axis equal
axis([-pi pi -pi pi])
end


function S=MakeSofR(R,M,n)
t=linspace(0,2*pi,n+1)';
t=t(1:end-1);
R=R*ones(M,1);
x0=2*pi*rand(M,1);
y0=2*pi*rand(M,1);
Sx=reshape(cos(t)*R'+ones(n,1)*x0',n*M,1);
Sy=reshape(sin(t)*R'+ones(n,1)*y0',n*M,1);
S=[Sx,Sy];
end

function S=MakeBarnettShape(R,M,n)
t=linspace(0,2*pi,n+1)';
t=t(1:end-1);
t=t*ones(1,M);
alpha=0.5*ones(n,1)*rand(1,M);
beta=2*pi*ones(n,1)*rand(1,M);
w=floor(7*ones(n,1)*rand(1,M)+2);
R=R*(1+alpha.*cos(w.*t+beta));
x=cos(t).*R;
y=sin(t).*R;
x0=2*pi*ones(n,1)*rand(1,M);
y0=2*pi*ones(n,1)*rand(1,M);
Sx=reshape(x+x0,n*M,1);
Sy=reshape(y+y0,n*M,1);
S=[Sx,Sy];
end