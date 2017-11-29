function u=Lap_SLP_Self_Apply(s,f,S)

M=length(s.len);
u.s=zeros(s.cumm(end),1);
ls=0;
for k=1:M
    ss.x=s.x(ls+1:ls+s.len(k));
    ss=quadr(ss);
    ff=f(ls+1:ls+s.len(k));
    if nargin==2
        zz=LapSLPselfmatrix(ss)*ff;
    else
        zz=S{k}*ff;
    end
    SLP=LapSLPmatrix(ss,ss,0);
    SLP(logical(eye(s.len(k))))=0;
    SLP=SLP*ff;
    u.s(ls+1:ls+s.len(k))=zz-SLP;
    ls=ls+s.len(k);
end