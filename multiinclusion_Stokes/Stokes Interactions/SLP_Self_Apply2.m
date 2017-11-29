function u=SLP_Self_Apply2(s,f,S)

M=length(s.len);
u.s=zeros(s.cumm(end),1);
ls=0;
for k=1:M
    ss.X=s.x([ls+s.len(k),ls+1:ls+s.len(k)]');
    ss.x=s.x(ls+1:ls+s.len(k));
    ss.w=s.w(ls+1:ls+s.len(k));
    ff=f(ls+1:ls+s.len(k));
    if nargin==2
        zz=SLPselfmatrix([real(ss.X),imag(ss.X)])*[real(ff);imag(ff)];
    else
        zz=S{k}*[real(ff);imag(ff)];
    end
    SLP=SLPmatrix(ss,ss,1);
    SLP(logical([eye(s.len(k)),eye(s.len(k));eye(s.len(k)),eye(s.len(k))]))=0;
    SLP=SLP*[real(ff);imag(ff)];
    zz=zz-SLP;
    u.s(ls+1:ls+s.len(k))=zz(1:end/2)+1i*zz(end/2+1:end);
    ls=ls+s.len(k);
end