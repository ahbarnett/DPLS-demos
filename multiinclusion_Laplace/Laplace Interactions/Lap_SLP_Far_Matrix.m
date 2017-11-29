function S=Lap_SLP_Far_Matrix(s,if_s,t,if_t)

M=length(s.len);

if if_s==1
    S.s=zeros(sum(s.len));
    S.sn=S.s;
    rs=0;
    for k=1:M
        ls=0;
        for j=1:M
            ss.x=s.x(rs+1:rs+s.len(k));
            ss.w=s.w(rs+1:rs+s.len(k));
            ss.nx=s.nx(rs+1:rs+s.len(k));
            ss.cur=s.cur(rs+1:rs+s.len(k));
            tt.x=s.x(ls+1:ls+s.len(j));
            tt.nx=s.nx(ls+1:ls+s.len(j));
   
            [SLP,SLPn]=LapSLPmatrix(tt,ss,0);
            if j==k
                SLP(logical(eye(s.len(k))))=0;
            end
            S.s(ls+1:ls+s.len(j),rs+1:rs+s.len(k))=SLP;
            S.sn(ls+1:ls+s.len(j),rs+1:rs+s.len(k))=SLPn;
            ls=ls+s.len(j);
        end
        rs=rs+s.len(k);
    end
end

if if_t==1
    S.t=zeros(length(t.x),sum(s.len));
    S.tx=S.t;
    S.ty=S.t;
    rs=0;
    for k=1:M
        ss.x=s.x(rs+1:rs+s.len(k));
        ss.w=s.w(rs+1:rs+s.len(k));
        SLP=LapSLPmatrix(t,ss,0);
        S.t(:,rs+1:rs+s.len(k))=SLP;
        
        t.nx=0*t.x+1;
        [~,SLPx]=LapSLPmatrix(t,ss,0);
        S.tx(:,rs+1:rs+s.len(k))=SLPx;
        
        t.nx=0*t.x+1i;
        [~,SLPy]=LapSLPmatrix(t,ss,0);
        S.ty(:,rs+1:rs+s.len(k))=SLPy;
        
        rs=rs+s.len(k);
    end
end