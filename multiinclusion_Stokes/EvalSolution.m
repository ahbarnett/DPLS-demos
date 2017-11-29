function u=EvalSolution(FileName,z)

load(FileName)
rng(RANDNUM)
warning off

uc.d = 2*pi; % unitcell, d=period in space

for k=1:MI
    I{k}=quadr(I{k});
end
mu = 0.7; % viscosity
uc.nei = 1; % how many nei copies either side (use eg 1e3 to test A w/o AP)
% n = 32; % pts per side
[x,w] = gauss(n); x = (1+x)/2; w = w'/2; % quadr on [0,1]
U.x=uc.d*x-uc.d/2+1i*uc.d/2; U.w=uc.d*w; U.nx=0*U.x-1i;
D=U; D.x=U.x-1i*uc.d;
H = uc.d; L.x = -uc.d/2+ 1i*H*x-uc.d/2*1i; L.nx = 0*L.x+1; L.w = H*w; % left side
R = L; R.x = L.x+uc.d; % right side
b.x = 1.4*uc.d*exp(2i*pi*(1:M)'/M);        % the proxy pts
b.w = 1+0*b.x; % unit quadr weights (dummy)

s.x=zeros(sum(N),1);
ls=0;
h=0;
for k=1:MI
    s.x(ls+1:ls+N(k))=I{k}.x;
    ls=ls+N(k);
    h0=ArcLength(real(I{k}.x([end,1:end])),imag(I{k}.x([end,1:end])),2*pi)/N(k);
    h=max(h0,h);
end
s.len=N;
s=Quad(s);

n=length(z.x);
u=zeros(2*n,1);
ls=0;
while ls<n
    lsmax=min(ls+2e5,n);
    zz.x=z.x(ls+1:lsmax);
    uu=At_operator(Eta,I,N,s,zz,mu,h);
    u([ls+1:lsmax,end/2+ls+1:end/2+lsmax])=u([ls+1:lsmax,end/2+ls+1:end/2+lsmax])+SLPmatrix(zz,b,mu)*Xi+uu;
    ls=lsmax;
end
end