function C=MakeC(s,L,R,T,B,U,nei)
M=length(L.x);
C=zeros(4*M,sum(s.len));
ls=0;
for k=1:s.M
    ss.x=s.x(ls+1:ls+s.len(k));
    ss=quadr(ss);
    N=s.len(k);
    CL = zeros(M,N); CLn = CL;  % fill C...
    for i=-nei:nei
        [Ci,Cni] = LapSLPmatrix(L,ss,nei*U.e1+i*U.e2); CL = CL+Ci; CLn=CLn+Cni;
    end
    CR = zeros(M,N); CRn = CR;
    for i=-nei:nei
        [Ci,Cni] = LapSLPmatrix(R,ss,-nei*U.e1+i*U.e2); CR = CR+Ci; CRn=CRn+Cni;
    end
    CT = zeros(M,N); CTn = CT;
    for i=-nei:nei
        [Ci,Cni] = LapSLPmatrix(T,ss,-nei*U.e2+i*U.e1); CT = CT+Ci; CTn=CTn+Cni;
    end
    CB = zeros(M,N); CBn = CB;
    for i=-nei:nei
        [Ci,Cni] = LapSLPmatrix(B,ss,nei*U.e2+i*U.e1); CB = CB+Ci; CBn=CBn+Cni;
    end
    C(:,ls+1:ls+s.len(k)) = [CR-CL; CRn-CLn; CT-CB; CTn-CBn];
    ls=ls+s.len(k);
end