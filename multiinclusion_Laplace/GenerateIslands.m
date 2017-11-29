function [I,N]=GenerateIslands(M,n)
% n = desired Nk

S=[];   % will be Nx2 list of all points of curves as built up
R0=10*sqrt(10)/6.9/sqrt(M);
rate=1.05;
R=R0*rate;
n0 = 1e3;    % Nk for building geom only
while size(S,1)/n0<M
    R=R/rate;
    MR=round(5/R^2);
    %S0=MakeSofR(R,MR,n);
    S0=MakeBarnettShape(R,MR,n0);
    S0=IslandIntersect(S0,S,n0);
    S=IslandIntersect(S,S0,n0);
    S0=RecursiveIntersect(S0,n0);
    S0=RecursiveIntersect(S0,n0);
    S=[S;S0];
end
S=S(1:M*n0,:);

ls=0;
I=cell(M,1);
for k=1:M
    I{k}.x=S(ls+1:ls+n0,:)*[1;1i]-pi-1i*pi;  % convert to C-# and translate
    area=Area(real(I{k}.x([1:end,1])),imag(I{k}.x([1:end,1])));  % true area
    xc=Area(real(I{k}.x([1:end,1])).^2/2,imag(I{k}.x([1:end,1])))/area;
    yc=Area(real(I{k}.x([1:end,1])).*imag(I{k}.x([1:end,1])),imag(I{k}.x([1:end,1])))/area;  % mystery - centroid of area assuming cont mass density?
    I{k}.x=0.97*(I{k}.x-xc-1i*yc)+xc+1i*yc;   % rescale about (xc,yc)
    % 0.97
    I{k}.x = perispecinterp(I{k}.x,n);    % resample to desired N_k, ahb
    ls=ls+n0;
end
N=n*ones(length(I),1);

for k=1:M
    for i=-1:1
        for j=-1:1
            plot(I{k}.x([1:end,1])+2*pi*i+1i*2*pi*j)
            hold on
        end
    end
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
%w=floor(7*ones(n,1)*rand(1,M)+2);
w=floor(4*ones(n,1)*rand(1,M)+2);
R=R*(1+alpha.*cos(w.*t+beta));
x=cos(t).*R;
y=sin(t).*R;
x0=2*pi*ones(n,1)*rand(1,M);
y0=2*pi*ones(n,1)*rand(1,M);
Sx=reshape(x+x0,n*M,1);
Sy=reshape(y+y0,n*M,1);
S=[Sx,Sy];
end