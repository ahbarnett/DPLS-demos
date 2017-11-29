function flux1=SolutionFlux(FileName)
% This function computes the flux between the left and right sides of the
% domain. 

load(FileName)
if ~exist('P','var')
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

proxyrep = @LapSLPmatrix;

U.e1 = 1; U.e2 = 1i; U.nei = 0; 
[x,w] = gauss(M); w=w/2;  % walls
L.x = (-U.e1 + U.e2*x)/2; L.nx = (-1i*U.e2)/abs(U.e2) + 0*L.x; L.w=w*abs(U.e2);
R = L; R.x = L.x + U.e1;
B.x = (-U.e2 + U.e1*x)/2; B.nx = (1i*U.e1)/abs(U.e1) + 0*B.x; B.w=w*abs(U.e1);
T = B; T.x = B.x + U.e2;
Rp = 1.4;
p.x = Rp * exp(1i*(1:P)'/P*2*pi); p = quadr(p); % proxy pts
U.L=L;
U.R=R;
U.T=T;
U.B=B;

w=[]; w.x = [L.x]; w.nx = [L.nx]; % flux through L wall
[~,Tw] = proxyrep(w, p); u = Tw*Psi;  % proxy contrib to L wall only
% this is crude, sums up u values on L-type and B-type walls
w=[]; w.nx = B.nx;
for i=-nei:nei   % top and bottom walls of 3x3
  w.x = B.x+2*U.e2; [~,DTw] = LapSLPmatrix(w, s, -U.e1*i);    % top
  u = u + (i+2)*DTw*Sig; % translate & integer factor  
  w.x = B.x-U.e2; [~,DTw] = LapSLPmatrix(w, s, -U.e1*i);    % bottom (flip nor)
  u = u - (i+2)*DTw*Sig; % translate & integer factor  
end
w=[]; w.nx = L.nx;
for j=-nei:nei     % right wall of 3x3
  w.x = L.x + j*U.e2;
  [~,DTw] = LapSLPmatrix(w, s, -U.e1*2); u = u + 3*DTw*Sig;
end
flux1 = L.w*u(1:M);  % assumes L and B have same quadr wei for now