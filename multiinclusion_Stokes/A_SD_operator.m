function Z1=A_SD_operator(X,I,N,s,S,t,h,mu,NPt,CSLP,CDLP)

MI=length(I);
X1=zeros(sum(N),1);
ls=0;
rs=0;
for k=1:MI
    X1(ls+1:ls+N(k))=X(rs+1:rs+N(k))+1i*X(rs+N(k)+1:rs+2*N(k));
    ls=ls+N(k);
    rs=rs+2*N(k);
end

Z=zeros(sum(N),1);
u_self=SLP_Self_Apply(s,X1,S);

if nargin==8
    NPt=FindNearPt(s,1,t,1,5*h);
end
if nargin<=9
    u_S_close=SLP_Close_Apply(NPt,s,t,X1);
    u_D_close=DLP_Close_Apply(NPt,s,t,X1);
else
    u_S_close=SLP_Close_Apply(NPt,s,t,X1,CSLP);
    u_D_close=DLP_Close_Apply(NPt,s,t,X1,CDLP);
end
u_S_far=SLP_Far_Apply(s,1,t,1,X1);
u_D_far=DLP_Far_Apply(s,1,t,1,X1);
Z=Z+u_self.s/mu;
Z=Z+u_S_close.s/mu+u_D_close.s+u_S_far.s/mu+u_D_far.s;

ls=0;
for k=1:8
Z=Z+u_S_close.t(ls+1:ls+sum(N))/mu+u_D_close.t(ls+1:ls+sum(N))+u_S_far.t(ls+1:ls+sum(N))/mu+u_D_far.t(ls+1:ls+sum(N));
ls=ls+sum(N);
end

Z1=zeros(2*sum(N),1);
ls=0;
rs=0;
for k=1:MI
    Z1(ls+1:ls+2*N(k))=[real(Z(rs+1:rs+N(k)));imag(Z(rs+1:rs+N(k)))];
    ls=ls+2*N(k);
    rs=rs+N(k);
end
Z1=1/2*X+Z1;