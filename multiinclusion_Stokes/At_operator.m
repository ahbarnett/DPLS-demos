function Z=At_operator(X,I,N,s,t,mu,h,NPt)

% s = source struct, t = target struct

MI=length(I);
X1=zeros(sum(N),1);
M=length(t.x);
ls=0;
rs=0;
for k=1:MI
    X1(ls+1:ls+N(k))=X(rs+1:rs+N(k))+1i*X(rs+N(k)+1:rs+2*N(k));
    ls=ls+N(k);
    rs=rs+2*N(k);
end

t0.x=[];
for i=-1:1
    for j=-1:1
        t0.x=[t0.x;t.x+2*pi*i+2*pi*1i*j];
    end
end
t.x=t0.x;

% S=SLP_Self_Matrix(s);
Z=zeros(M,1);
% u_self=SLP_Self_Apply(s,X1,S);
%u_self=SLP_Self_Apply(s,X1);

if nargin==7
    NPt=FindNearPt(s,0,t,1,5*h);
end
u_S_close=SLP_Close_Apply(NPt,s,t,X1);
u_D_close=DLP_Close_Apply(NPt,s,t,X1);
u_S_far=SLP_Far_Apply(s,0,t,1,X1);
u_D_far=DLP_Far_Apply(s,0,t,1,X1);
% Z=Z+u_self.s/mu;
% Z=Z+u_S_close.s/mu+u_D_close.s+u_S_far.s/mu+u_D_far.s;

ls=0;
for k=1:9
Z=Z+u_S_close.t(ls+1:ls+M)/mu+u_D_close.t(ls+1:ls+M)+u_S_far.t(ls+1:ls+M)/mu+u_D_far.t(ls+1:ls+M);
ls=ls+M;
end

Z=[real(Z);imag(Z)];

% A = zeros(2*sum(N)); % A's jump relation part.
% ls=0;
% for k=1:MI
%     for j=-1:1
%         for i=-1:1
%             I0=I{k};
%             I0.x=I{k}.x+uc.d*1i*j+uc.d*i;
%             A(ls+1:ls+2*N(k),ls+1:ls+2*N(k))=A(ls+1:ls+2*N(k),ls+1:ls+2*N(k))+SLPmatrix_scf(I{k}.x,I0,mu);
%             for kk=[1:k-1,k+1:MI]
%                 I1.x=I{kk}.x+uc.d*1i*j+uc.d*i;
%                 I1=quadr(I1);
%                 lk=2*sum(N(1:kk-1));
%                 A(ls+1:ls+2*N(k),lk+1:lk+2*N(kk))=...
%                     A(ls+1:ls+2*N(k),lk+1:lk+2*N(kk))+...
%                     SLPmatrix_scf(I{k}.x,I1,mu);
%             end
%         end
%     end
%    ls=ls+2*N(k);
% end
% Z=A*X;
% max(abs(Z-Z1))