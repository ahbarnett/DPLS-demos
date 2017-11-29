function [u,ux,uy]=At_operator(t,s,f,U,nei)

m=length(t.x);
tt.x=zeros(9*m,1);
ls=0;
for i=-nei:nei
    for j=-nei:nei
            tt.x(ls+1:ls+m)=t.x+U.e1*i+U.e2*j;
            ls=ls+m;
    end
end


NPt=FindNearPt(s,0,tt,1,5*max(s.h));
u_close=Lap_SLP_Close_Apply(NPt,s,tt,f);
u_far=Lap_SLP_Far_Apply(s,0,tt,1,f);

u=zeros(m,1);
ux=u;
uy=u;
ls=0;
for j=1:9
   u=u+u_close.t(ls+1:ls+m)+u_far.t(ls+1:ls+m);
   ux=ux+u_close.tx(ls+1:ls+m)+u_far.tx(ls+1:ls+m);
   uy=uy+u_close.ty(ls+1:ls+m)+u_far.ty(ls+1:ls+m);
   ls=ls+m;
end