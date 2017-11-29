function u=A_operator(s,t,f,NPt)

n=sum(s.len);
if nargin==3
    NPt=FindNearPt(s,1,t,1,5*max(s.h));
end
u_close=Lap_SLP_Close_Apply(NPt,s,t,f);
u_far=Lap_SLP_Far_Apply(s,1,t,1,f);

u=u_close.sn+u_far.sn;
ls=0;
for j=1:8
   u=u+(u_close.tx(ls+1:ls+n)+u_far.tx(ls+1:ls+n)).*real(s.nx)+(u_close.ty(ls+1:ls+n)+u_far.ty(ls+1:ls+n)).*imag(s.nx);
   ls=ls+n;
end