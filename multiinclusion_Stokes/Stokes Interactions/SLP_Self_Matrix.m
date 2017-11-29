function S=SLP_Self_Matrix(s)

M=length(s.len);
S=cell(1,M);
ls=0;
for k=1:M
   S{k}=SLPselfmatrix([real(s.x([ls+s.len(k),ls+1:ls+s.len(k)]')),imag(s.x([ls+s.len(k),ls+1:ls+s.len(k)]'))]);
   ls=ls+s.len(k);
end