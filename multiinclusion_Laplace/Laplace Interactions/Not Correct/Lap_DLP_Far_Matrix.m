function D=Lap_DLP_Far_Matrix(s,if_s,t,if_t)

M=length(s.len);

if if_s==1
    D.s=zeros(sum(s.len));
    rs=0;
    for k=1:M
        ls=0;
        for j=1:M
            ss.x=s.x(rs+1:rs+s.len(k));
            ss.nx=s.nx(rs+1:rs+s.len(k));
            ss.cur=s.cur(rs+1:rs+s.len(k));
            ss.w=s.w(rs+1:rs+s.len(k));
            tt.x=s.x(ls+1:ls+s.len(j));
            DLP=LapDLPmatrix(tt,ss,0);
            D.s(ls+1:ls+s.len(j),rs+1:rs+s.len(k))=DLP;
            ls=ls+s.len(j);
        end
        rs=rs+s.len(k);
    end
end

if if_t==1
    D.t=zeros(length(t.x),sum(s.len));
    rs=0;
    for k=1:M
        ss.x=s.x(rs+1:rs+s.len(k));
        ss.nx=s.nx(rs+1:rs+s.len(k));
        ss.cur=s.cur(rs+1:rs+s.len(k));
        ss.w=s.w(rs+1:rs+s.len(k));
        DLP=LapDLPmatrix(t,ss,0);
        D.t(:,rs+1:rs+s.len(k))=DLP;
        rs=rs+s.len(k);
    end
end