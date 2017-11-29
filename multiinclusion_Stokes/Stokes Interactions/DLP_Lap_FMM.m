function [u,ux,uy]=DLP_Lap_FMM(s,if_s,t,if_t,sig,iprec)

if nargin==5
    iprec=4;
end

nsource=length(s.x);
source=[real(s.x),imag(s.x)]';
dipstr = (-1i*s.xpn.*sig).';
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

U=zfmm2dpart(iprec,nsource,source,dipstr,ifpot,ifgrad,ifhess,ntarget,target,ifpottarg,ifgradtarg,ifhesstarg);

if if_s==1
    u.s=real(U.pot.');
    ux.s=zeros(nsource,1);
    uy.s=ux.s;
    ux.s(1:2:end)=real(U.grad(1,1:end/2).');
    ux.s(2:2:end)=real(U.grad(2,1:end/2).');
    uy.s(1:2:end)=-imag(U.grad(1,1:end/2).');
    uy.s(2:2:end)=-imag(U.grad(2,1:end/2).');
end

if if_t==1
    u.t=real(U.pottarg.');
    ux.t=zeros(ntarget,1);
    uy.t=ux.t;
    ux.t(1:2:end)=real(U.gradtarg(1,1:ceil(end/2)).');
    ux.t(2:2:end)=real(U.gradtarg(2,1:floor(end/2)).');
    uy.t(1:2:end)=-imag(U.gradtarg(1,1:ceil(end/2)).');
    uy.t(2:2:end)=-imag(U.gradtarg(2,1:floor(end/2)).');
end