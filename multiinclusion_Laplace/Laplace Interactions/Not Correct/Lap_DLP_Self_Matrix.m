function D=Lap_DLP_Self_Matrix(s)

M=length(s.len);
D=cell(1,M);
ls=0;
for k=1:M
   ss.x=s.x(ls+1:ls+s.len(k));
   ss.nx=s.nx(ls+1:ls+s.len(k));
   ss.cur=s.cur(ls+1:ls+s.len(k));
   ss.w=s.w(ls+1:ls+s.len(k));
   D{k}=LapDLPmatrix(ss,ss,0);
   ls=ls+s.len(k);
end