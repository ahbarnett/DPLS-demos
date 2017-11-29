function s=Quad(s)

s.M=length(s.len);
s.cumm=s.len;
for k=2:s.M
   s.cumm(k)=s.cumm(k)+s.cumm(k-1);
end
s.t=0*s.x(:);
s.xp=s.t;
s.xpp=s.t;
s.sp=s.t;
s.tang=s.t;
s.nx=s.t;
s.cur=s.t;
s.w=s.t;
s.cw=s.t;
s.xpn=s.t;
s.spn=s.t;
s.h=0*s.len;
ls=0;
for k=1:s.M
    n=s.len(k);
    ss.x=s.x(ls+1:ls+n);
    ss=quadr(ss);
    s.xp(ls+1:ls+n)=ss.xp;
    s.xpp(ls+1:ls+n)=ss.xpp;
    s.t(ls+1:ls+n)=ss.t;
    s.sp(ls+1:ls+n)=ss.sp;
    s.tang(ls+1:ls+n)=ss.tang;
    s.nx(ls+1:ls+n)=ss.nx;
    s.cur(ls+1:ls+n)=ss.cur;
    s.w(ls+1:ls+n)=ss.w;
    s.cw(ls+1:ls+n)=ss.cw;
    s.xpn(ls+1:ls+n)=ss.xp/n;
    s.spn(ls+1:ls+n)=ss.sp/n;
    s.h(k)=max(2*pi*ss.sp/n);
    ls=ls+n;
end