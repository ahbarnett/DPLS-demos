function LaplaceSolver(FileName)
% Solve 2D Neumann Laplace doubly-periodic BVP, via proxy periodizing

% Load parameters and inclusion positions.
load(FileName)
if ~exist('P','var')
    P=M;
    M=n;
    clear n
end
nei = 1;      % 0 for no direct neighbors; 1 for preferred 3x3 scheme
jumps = [1 0]; %[0 1]; % potential jumps across R-L and T-B

% Computes inclusion properties
s.x=zeros(sum(N),1);
s.len=N;
ls=0;
for k=1:length(N)
    s.x(ls+1:ls+s.len(k))=I{k}.x/2/pi;
    ls=ls+s.len(k);
end
s=Quad(s);

tic
% Constructs sides of domain and proxy points
U.e1 = 1; U.e2 = 1i; U.nei = 0;  % unit cell
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

% Builds matrices/operators
t.x=zeros(8*sum(s.len),1);
t.nx=t.x;
ls=0;
for i=-nei:nei
    for j=-nei:nei
        if i~=0||j~=0
            t.x(ls+1:ls+sum(s.len))=s.x+U.e1*i+U.e2*j;
            t.nx(ls+1:ls+sum(s.len))=s.nx;
            ls=ls+sum(s.len);
        end
    end
end

NPt=FindNearPt(s,1,t,1,5*max(s.h)); % Determines which points are within 5h of an island
Aop=@(x)-x/2+A_operator(s,t,x,NPt); % A operator
[~,Bm] = LapSLPmatrix(s,p);    % B operator
C=MakeC(s,L,R,T,B,U,nei); % C operator
[QL,QLn] = LapSLPmatrix(L,p); [QR,QRn] = LapSLPmatrix(R,p);
[QB,QBn] = LapSLPmatrix(B,p); [QT,QTn] = LapSLPmatrix(T,p);
Q = [QR-QL; QRn-QLn; QT-QB; QTn-QBn]; % Q operator

% Constructs right-hand side vector
rhs = [0*s.x; jumps(1)+0*L.x; 0*L.x; jumps(2)+0*B.x; 0*B.x]; % box driving    
f=rhs(1:sum(s.len)); g=rhs(sum(s.len)+1:end);

% Expands the column space of Q
H=0*s.x+1;
R=ones(P,1)/P;
Qv=Q+(C*H)*R'; % Qv is Q with an expanded column space

% Applying Qv^+ to g and C
lso.RECT = true;
Qvdagg=linsolve(Qv,g,lso);
QvdagC=linsolve(Qv,C,lso);

PreT=toc;
tic

% Uses gmres to solve (A-B(Qv^+C))Sig=f-B(Qv^+g) 
[Sig,~,~,ITER] = gmres(@(x) Aop(x)-Bm*(QvdagC*x), f-Bm*Qvdagg,[],ACC,MAXIT);
GMRESTime=toc;
tic
% Computes Psi using Psi=Qv^+g-(Qv^+C)Sig
Psi = Qvdagg - QvdagC*Sig; % get remaining unknowns
PsiTime=toc;
Time=PreT+GMRESTime+PsiTime;

% Computes residual
Res=max(abs([Aop(Sig)+Bm*Psi;C*Sig+Q*Psi]-rhs));
ITER=ITER(2);

% Displays times
disp(['Elapsed time: ',num2str(Time),' seconds'])
disp(['GMRES iterations: ',num2str(ITER)])
disp(['Residual error: ',num2str(Res)])

% Saves variables
save(FileName,'Sig','Psi','Res','ITER','Time','GMRESTime','-append')
end
