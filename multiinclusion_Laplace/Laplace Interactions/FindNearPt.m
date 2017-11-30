function NPt=FindNearPt(s,if_s,t,if_t,d)

% hacked Barnett version which replaces Stats Toolbox kdtree
% with an open-source version;
% but I don't really understand Gary's code here. 11/30/17

Y=[real(s.x),imag(s.x)];
M=length(s.len);

if if_s==1       % self search over s
    NPt.s=cell(1,M);
    tree = kdtree_build(Y);   % Y is Mx2
    ls=0;
    for k=1:M
      idy = [];
      for m=ls+1:ls+s.len(k)   % can only query one at a time
        idy = [idy; kdtree_ball_query(tree, Y(m,:),d)];  % ahb
      end
      NPt.s{k}=unique(idy);
      NPt.s{k}=NPt.s{k}(NPt.s{k}<ls+1|NPt.s{k}>ls+s.len(k));
      ls=ls+s.len(k);
    end
end

if if_t==1      % tree on targs, search each source
    NPt.t=cell(1,M);
    X=[real(t.x),imag(t.x)];
    tree = kdtree_build(X);
    ls=0;
    for k=1:M
        idx = [];
        for m=ls+1:ls+s.len(k)   % can only query one at a time
          idx = [idx; kdtree_ball_query(tree, Y(m,:),d)];  % ahb
        end
        NPt.t{k}=unique(idx);
        ls=ls+s.len(k);
    end
end
