function flux1=SolutionFlux(FileName)
% This function computes the flux between the left and right sides of the
% domain.

load(FileName)
nei = 1;

uc.d = 2*pi;
mu = 0.7; % viscosity
uc.nei = 1; % how many nei copies either side (use eg 1e3 to test A w/o AP)
[x,w] = gauss(n); x = (1+x)/2; w = w'/2; % quadr on [0,1]
U.x=uc.d*x-uc.d/2+1i*uc.d/2; U.w=uc.d*w; U.nx=0*U.x-1i;
D=U; D.x=U.x-1i*uc.d;
H = uc.d;
D.w=H*w;
L.x = -uc.d/2+ 1i*H*x-uc.d/2*1i; L.nx = 0*L.x+1; L.w = H*w; % left side
R = L; R.x = L.x+uc.d; % right side
b.x = 1.4*uc.d*exp(2i*pi*(1:M)'/M);        % the proxy pts
b.w = 1+0*b.x; % unit quadr weights (dummy)

w=[]; w.x=L.x; w.nx=L.nx; % flux through L wall
u=SLPmatrix(w,b,mu)*Xi;
u0=u(1:end/2);
u=0;
w=[]; w.nx = D.nx;
ls=0;
for k=1:MI
    I0=I{k};
    I0=quadr(I0);
    for i=-nei:nei   % top and bottom walls of 3x3
        w.x = D.x+4*pi*1i; 
        u=u+(i+2)*(SLPmatrix(w,I0,mu,-2*pi*i)+DLPmatrix(w,I0,mu,-2*pi*i))*Eta(ls+1:ls+2*N(k));    % top
        w.x = D.x-2*pi*1i; 
        u=u-(i+2)*(SLPmatrix(w,I0,mu,-2*pi*i)+DLPmatrix(w,I0,mu,-2*pi*i))*Eta(ls+1:ls+2*N(k));  
    end
    ls=ls+2*N(k);
end
u0=u0+u(end/2+1:end);
u=0;
w=[]; w.nx = L.nx;
ls=0;
for k=1:MI
    I0=I{k};
    I0=quadr(I0);
    for j=-nei:nei     % right wall of 3x3
        w.x = L.x + j*2*pi*1i;
        u=u+3*(SLPmatrix(w,I0,mu,-4*pi)+DLPmatrix(w,I0,mu,-4*pi))*Eta(ls+1:ls+2*N(k));
    end
    ls=ls+2*N(k);
end
u0=u0+u(1:end/2);
flux1 = L.w'*u0;  % assumes L and B have same quadr wei for now