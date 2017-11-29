function D=DLP_Far_Matrix(s,if_s,t,if_t)

M=length(s.len);

if if_s==1
    D.s=zeros(2*sum(s.len));
    rs=0;
    for k=1:M
        ls=0;
        for j=1:M
            ss.x=s.x(rs+1:rs+s.len(k));
            ss.nx=s.nx(rs+1:rs+s.len(k));
            ss.cur=s.cur(rs+1:rs+s.len(k));
            ss.w=s.w(rs+1:rs+s.len(k));
            tt.x=s.x(ls+1:ls+s.len(j));
            DLP=DLPmatrix(tt,ss,1);
            D.s(ls+1:ls+s.len(j),rs+1:rs+s.len(k))=DLP(1:end/2,1:end/2);
            D.s(end/2+ls+1:end/2+ls+s.len(j),rs+1:rs+s.len(k))=DLP(end/2+1:end,1:end/2);
            D.s(ls+1:ls+s.len(j),end/2+rs+1:end/2+rs+s.len(k))=DLP(1:end/2,end/2+1:end);
            D.s(end/2+ls+1:end/2+ls+s.len(j),end/2+rs+1:end/2+rs+s.len(k))=DLP(end/2+1:end,end/2+1:end);
            ls=ls+s.len(j);
        end
        rs=rs+s.len(k);
    end
end

if if_t==1
    D.t=zeros(2*length(t.x),2*sum(s.len));
    rs=0;
    for k=1:M
        ss.x=s.x(rs+1:rs+s.len(k));
        ss.nx=s.nx(rs+1:rs+s.len(k));
        ss.cur=s.cur(rs+1:rs+s.len(k));
        ss.w=s.w(rs+1:rs+s.len(k));
        DLP=DLPmatrix(t,ss,1);
        D.t(1:end/2,rs+1:rs+s.len(k))=DLP(1:end/2,1:end/2);
        D.t(end/2+1:end,rs+1:rs+s.len(k))=DLP(end/2+1:end,1:end/2);
        D.t(1:end/2,end/2+rs+1:end/2+rs+s.len(k))=DLP(1:end/2,end/2+1:end);
        D.t(end/2+1:end,end/2+rs+1:end/2+rs+s.len(k))=DLP(end/2+1:end,end/2+1:end);
        rs=rs+s.len(k);
    end
end