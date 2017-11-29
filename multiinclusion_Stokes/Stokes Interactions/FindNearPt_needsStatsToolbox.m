function NPt=FindNearPt(s,if_s,t,if_t,d)

Y=[real(s.x),imag(s.x)];
M=length(s.len);

if if_s==1
    NPt.s=cell(1,M);
    ny=createns(Y,'nsmethod','kdtree');
    idy=rangesearch(ny,Y,d);
    idy=idy';
    ls=0;
    for k=1:M
        NPt.s{k}=idy(ls+1:ls+s.len(k));
        NPt.s{k}=unique(cell2mat(NPt.s{k}))';
        NPt.s{k}=NPt.s{k}(NPt.s{k}<ls+1|NPt.s{k}>ls+s.len(k));
        ls=ls+s.len(k);
    end
end

if if_t==1
    NPt.t=cell(1,M);
    X=[real(t.x),imag(t.x)];
    nx=createns(X,'nsmethod','kdtree');
    idx=rangesearch(nx,Y,d);
    idx=idx';
    ls=0;
    for k=1:M
        NPt.t{k}=idx(ls+1:ls+s.len(k));
        NPt.t{k}=unique(cell2mat(NPt.t{k}))';
        ls=ls+s.len(k);
    end
end