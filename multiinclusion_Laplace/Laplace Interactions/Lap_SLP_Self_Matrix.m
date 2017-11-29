function [S,Sn]=Lap_SLP_Self_Matrix(s)

M=length(s.len);
S=cell(1,M);
Sn=S;
ls=0;
for k=1:M
    ss.x=s.x(ls+1:ls+s.len(k));
    ss=quadr(ss);
    S{k}=LapSLPselfmatrix(ss);
    if nargout==2
        Sn{k}=LapSLPmatrix(ss,ss,0);
    end
    ls=ls+s.len(k);
end