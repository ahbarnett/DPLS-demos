function StokesSolver(FileName,MEM)
% Doubly periodic Stokes flow with "island" geom using circle of SLP proxy
% sources. Gary Marple using code from Alex Barnett 6/5/15. If plent of
% memory is available, set MEM=1 to create and store self-interaction
% matrices to speed up GMRES iterations.

load(FileName)
rng(RANDNUM)
warning off

tic
uc.d = 2*pi; % unitcell, d=period in space

% Computes quadrature information for inclusions
for k=1:MI
    I{k}=quadr(I{k});
end
mu = 0.7; % viscosity

% Constructs sides of domain and proxy points
uc.nei = 1; % how many nei copies either side (use eg 1e3 to test A w/o AP)
[x,w] = gauss(n); x = (1+x)/2; w = w'/2;
U.x=uc.d*x-uc.d/2+1i*uc.d/2; U.w=uc.d*w; U.nx=0*U.x-1i;
D=U; D.x=U.x-1i*uc.d;
H = uc.d; L.x = -uc.d/2+ 1i*H*x-uc.d/2*1i; L.nx = 0*L.x+1; L.w = H*w; % left side
R = L; R.x = L.x+uc.d; % right side
b.x = 1.4*uc.d*exp(2i*pi*(1:M)'/M);        % the proxy pts
b.w = 1+0*b.x; % unit quadr weights (dummy)

% Forming the f vector
if 1
    f=zeros(2*sum(N),1); % No-slip BCs
else
    f=[];
    for k=1:MI
        f=[f;ue(I{k}.x)]; % Driving: I bdry vels, each stacked as [u_1;u_2]
    end
end

% Fomring the g vector
alpha=1; % traction jump in x-direction
beta=0; % traction jump in y-direction
g=[zeros(2*n,1);alpha*ones(n,1);zeros(4*n,1);beta*ones(n,1)];

% Construct self interactions matrices
tic
s.x=zeros(sum(N),1);
CSLP=cell(1,length(N));
CDLP=cell(1,length(N));
ls=0;
h=0;
for k=1:MI
    s.x(ls+1:ls+N(k))=I{k}.x;
    I{k}=quadr(I{k});
    if nargin==2&&MEM==1
        CSLP{k}=CSLPselfmatrix(I{k},'e');
        CDLP{k}=CDLPmatrix(I{k});
    end
    ls=ls+N(k);
    h0=ArcLength(real(I{k}.x([end,1:end])),imag(I{k}.x([end,1:end])),2*pi)/N(k);
    h=max(h0,h);
end
s.len=N;
s=Quad(s);
t.x=zeros(8*sum(N),1);
ls=0;
for j=-1:1
    for i=-1:1
        if i~=0||j~=0
            for k=1:MI
                t.x(ls+1:ls+s.len(k))=I{k}.x+i*uc.d+j*1i*uc.d;
                ls=ls+s.len(k);
            end
        end
    end
end

% Computes self matrices
SM=SLP_Self_Matrix(s);
ls=0;
for k=1:MI
    ss.x=s.x(ls+1:ls+s.len(k));
    ss.w=s.w(ls+1:ls+s.len(k));
    SLP=SLPmatrix(ss,ss,1);
    SLP(logical([eye(s.len(k)),eye(s.len(k));eye(s.len(k)),eye(s.len(k))]))=0;
    SM{k}=SM{k}-SLP;
    ls=ls+s.len(k);
end

% Determines which points are within 5h of an inclusion
NPt=FindNearPt(s,1,t,1,5*h);
if nargin==2&&MEM==1
    Aop=@(x)A_SD_operator(x,I,N,s,SM,t,h,mu,NPt,CSLP,CDLP);
else
    Aop=@(x)A_SD_operator(x,I,N,s,SM,t,h,mu,NPt);
end

% single layer (monopoles) on proxy points:
B=zeros(2*sum(N),2*M);
ls=0;
for k=1:MI
    B(ls+1:ls+2*N(k),:)=SLPmatrix(I{k},b,mu);
    ls=ls+2*N(k);
end

% Computes C
C=zeros(8*n,2*sum(N));
ls=0;
for k=1:MI
    for j=-1:1
        a=uc.nei*uc.d;
        I0=I{k}; I0.x=I{k}.x+uc.d*1i*j; IU=I{k}; ID=I{k};
        IU.x=I{k}.x+uc.d*j-uc.d*1i; ID.x=I{k}.x-uc.d*j+uc.d*1i;
        [S_RI,S_RIn] = SLPmatrix(R,I0,mu,-a); [S_LI,S_LIn] = SLPmatrix(L,I0,mu,a);
        [S_UI,S_UIn] = SLPmatrix(U,IU,mu); [S_DI,S_DIn] = SLPmatrix(D,ID,mu);
        [D_RI,D_RIn] = DLPmatrix(R,I0,mu,-a); [D_LI,D_LIn] = DLPmatrix(L,I0,mu,a);
        [D_UI,D_UIn] = DLPmatrix(U,IU,mu); [D_DI,D_DIn] = DLPmatrix(D,ID,mu);
        C(:,ls+1:ls+2*N(k)) = C(:,ls+1:ls+2*N(k))+...
            [S_RI-S_LI; S_RIn-S_LIn;S_DI-S_UI;S_DIn-S_UIn]+...
            [D_RI-D_LI; D_RIn-D_LIn;D_DI-D_UI;D_DIn-D_UIn];
        % maps cancelled densities to discrepancy
    end
    ls=ls+2*N(k);
end

% computes Q
[Rb,Rbn] = SLPmatrix(R,b,mu); [Lb,Lbn] = SLPmatrix(L,b,mu);
[Ub,Ubn] = SLPmatrix(U,b,mu); [Db,Dbn] = SLPmatrix(D,b,mu);
Q = [Rb-Lb; Rbn-Lbn;Db-Ub; Dbn-Ubn]; % maps periodizing dofs to discrepancy

% Expands the columns space of Q
H=zeros(2*sum(N),3);
ls=0;
rs=0;
for k=1:MI
H(ls+1:ls+2*N(k),:)=[ones(N(k),1),zeros(N(k),1),real(s.nx(rs+1:rs+N(k)));zeros(N(k),1),ones(N(k),1),imag(s.nx(rs+1:rs+N(k)))];
ls=ls+2*N(k);
rs=rs+N(k);
end
R=1/M*[ones(M,1),zeros(M,1),cos(2*pi*(1:M)'/M);zeros(M,1),ones(M,1),sin(2*pi*(1:M)'/M)];
Qv=Q+(C*H)*R';
Bw=B+([Aop(H(:,1)),Aop(H(:,2))])*R(:,1:2)';

% Computes Qv^+C and Qv^+g
lso.RECT = true;
QvdagC=linsolve(Qv,C,lso);
Qvdagg=linsolve(Qv,g,lso);

% Forms the right-hand side vector and the periodized A operator
rhs=f-Bw*Qvdagg;
Ap=@(x)Aop(x)-Bw*(QvdagC*x);

PreT=toc;
tic
[eta,~,~,ITER] = gmres(Ap,rhs,[],ACC,MAXIT); % Solves for eta using GMRES
GMRESTime=toc;
tic
xi = Qvdagg - QvdagC*eta; % Solves for xi


PsiTime=toc;
Time=PreT+GMRESTime+PsiTime;
ITER=ITER(2);
Eta=eta+H*(R'*xi);
Xi=xi;
Res=max(abs([Aop(Eta)+B*Xi-f;C*Eta+Q*Xi-g])); % Computes residual

% Displays timeing, gmres iterations, and the residual
disp(['Elapsed time: ',num2str(Time),' seconds'])
disp(['GMRES iterations: ',num2str(ITER)])
disp(['Residual error: ',num2str(Res)])

% Saves results
save(FileName,'Eta','Xi','Res','ITER','Time','GMRESTime','-append')
