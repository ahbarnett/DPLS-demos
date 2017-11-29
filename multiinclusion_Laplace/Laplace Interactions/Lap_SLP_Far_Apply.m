function u=Lap_SLP_Far_Apply(s,if_s,t,if_t,sig,iprec)

if nargin==5
    iprec=4;
end

nsource=length(s.x);
source=[real(s.x),imag(s.x)]';
ifcharge=1;
charge=-((s.spn).*sig)';
ifdipole=0;
dipstr = 0;
dipvec = zeros(2,nsource);
ifhess = 0;
ifhesstarg = 0;

if if_s==1
    ifpot=1;
    ifgrad=1;
else
    ifpot=0;
    ifgrad=0;
end

if if_t==1
    ntarget=length(t.x);
    target=[real(t.x),imag(t.x)]';
    ifpottarg=1;
    ifgradtarg=1;
else
    ntarget=0;
    target=zeros(2,1);
    ifpottarg=0;
    ifgradtarg=0;
end

U=rfmm2dpart(iprec,nsource,source,ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,ifgrad,ifhess,ntarget,target,ifpottarg,ifgradtarg,ifhesstarg);

if if_s==1
    u.s=U.pot.';
    u.sn=U.grad(1,:).'.*real(s.nx)+U.grad(2,:).'.*imag(s.nx);
    u.sn=u.sn-s.cur.*sig.*s.w/4/pi;
end

if if_t==1
    u.t=U.pottarg.';
    u.tx=U.gradtarg(1,:).';
    u.ty=U.gradtarg(2,:).';
end