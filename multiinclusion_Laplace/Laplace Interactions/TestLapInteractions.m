function TestLapInteractions

rng(1)
addpath(genpath('..'))
n=128;
alpha=linspace(0,2*pi,n+1)';
alpha=alpha(2:end);
s.x=cos(alpha)+1i*(sin(alpha)+0.5*cos(2*alpha));
s.len(1)=n;
%s.x=[s.x;cos(alpha)+1i*exp(sin(alpha))+0.23*1i];
%s.len(2)=n;
s=Quad(s);
f=[exp(sin(alpha));exp(cos(alpha))];
f=f(1:end/2);
t.x=(4*pi*rand(40,2)-2*pi)*[1;1i];
plot(real(s.x),imag(s.x))
hold on
scatter(real(t.x),imag(t.x),'r')
t.x=s.x+1e-10*s.nx;
t.nx=s.nx;
NPt=FindNearPt(s,1,t,1,100);
scatter(real(s.x(NPt.s{1})),imag(s.x(NPt.s{1})),'b');
%U=Lap_SLP_Close_Apply(NPt,s,t,f);
hold on
scatter(real(t.x),imag(t.x),'g')
[u,ux,uy]=lapSevalclose(t.x,s,f,'e');

U_Far=Lap_SLP_Far_Apply(s,1,t,1,f);
U_Close=Lap_SLP_Close_Apply(NPt,s,t,f);
S=Lap_SLP_Self_Matrix(s);
U_Self=Lap_SLP_Self_Apply(s,f,S);
SLP=LapSLPmatrix(s,s,0);
SLP(diag(true(size(SLP,1),1)))=0;
U_Self.s-(u-SLP*f)

U_Far.ty+U_Close.ty-uy;

SM=Lap_SLP_Far_Matrix(s,1,t,1);
[SM.sn*f-S.sn]




% U.sn
% [U.tx,U.ty]
% pause(100)
% size(s.x)
% size(f)
% size(t.x)
% 
% t.x=s.x;
% t.nx=s.nx*1i
% [A,An]=LapSLPmatrix(t,s,1);
% figure
% surf(An)
% pause(100)
% [~,ux,uy]=lapDevalclose(t.x,s,f,'e');
% [Uy.t,uy]
% V=Lap_DLP_Close_Apply(NPt,s,t,f);
% 
% 
% %S=Lap_SLP_Self_Matrix(s);
% %U_self=Lap_DLP_Self_Apply(s,f,S);
% U_far_M=Lap_SLP_Far_Matrix(s,1,t,1);
% U_far=Lap_SLP_Far_Apply(s,1,t,1,f);
% U_close=Lap_DLP_Close_Apply(NPt,s,t,f);
% ls=0;
% U_exact=zeros(length(t.x),1);
% for k=1:length(s.len)
% ss.x=s.x(ls+1:ls+s.len(k));
% ss=quadr(ss);
% %tt.x=s.x+1e2*eps*s.nx;
% U_exact=U_exact+lapDevalclose(t.x,ss,f(ls+1:ls+s.len(k)),'e'); 
% ls=ls+s.len(k);
% end
% 
% U_far_M.t*f-U_far.t
% U=U_far.t+U_close.t;
% %plot(U,'b')
% 
% % 
% k=1;
% ls=0;
% [~,ux,uy]=lapSevalclose(t.x,ss,f(ls+1:ls+s.len(k)),'e');
% [~,uux,uuy]=lapDevalclose(t.x,ss,f(ls+1:ls+s.len(k)),'e');
% t.nx=0*t.x+1;
% [~,Anx]=LapSLPmatrix(t,ss);
% [ux-Anx*f(ls+1:ls+s.len(k))]
% [~,AAnx]=LapDLPmatrix(t,ss);
% uux-AAnx*f(ls+1:ls+s.len(k))