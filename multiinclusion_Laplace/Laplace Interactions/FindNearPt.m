function NPt=FindNearPt(s,if_s,t,if_t,d)

% Barnett version which replaces Stats Toolbox kdtree with an open-source
% version.

Y=[real(s.x),imag(s.x)];
M=length(s.len);
addpath ~/numerics/kdtree/toolbox

if if_s==1
    NPt.s=cell(1,M);
    tree = kdtree_build(Y);   % Y is Mx2
  %  ny=createns(Y,'nsmethod','kdtree');
 %   idy=rangesearch(ny,Y,d);
  %  idy=idy';
    ls=0;
    for k=1:M
        idy = kdtree_ball_query(tree, Y(k,:),d);  % ahb
        NPt.s{k}=idy; %(ls+1:ls+s.len(k));
        NPt.s{k}=unique(NPt.s{k})';
        NPt.s{k}=NPt.s{k}(NPt.s{k}<ls+1|NPt.s{k}>ls+s.len(k));
        ls=ls+s.len(k);
    end
end

if if_t==1
    NPt.t=cell(1,M);
    X=[real(t.x),imag(t.x)];
    tree = kdtree_build(X);
   % nx=createns(X,'nsmethod','kdtree');
   % idx=rangesearch(nx,Y,d);
   % idx=idx';
    ls=0;
    for k=1:M
        idx = kdtree_ball_query(tree, Y(k,:),d);  % ahb
        NPt.t{k}=idx; %(ls+1:ls+s.len(k));
        NPt.t{k}=unique(NPt.t{k})';
        ls=ls+s.len(k);
    end
end