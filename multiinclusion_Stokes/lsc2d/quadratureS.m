function quadr = quadratureS(m, q, L)
%QUADRATURES gives nodes and weights to integrate logarithmic sigularity.
%
ker = 2; l1 = 0; l2 = L/2;
[v, u, a] = give_the_weights(q, ker);
[x, w] = regular_weights(a);

n = m-2*a+1;
h = 1/m;

evalpots1 = v*h;
evalpots3 = 1-x*h;

ys = l1 + (l2 - l1)*[evalpots1;evalpots3];
wt = (l2-l1)*h*[u; w];

wt = [wt; flipud(wt)]; ys = [ys; L-flipud(ys)];
A = sinterpS(2*m, ys, L); h = L/2/m; yt = [a*h:h:(m - a)*h]';
wt = [wt; h*ones(2*length(yt),1)]/4/pi;
lyt = 2*length(yt); B = sparse(lyt, 2*m); pos = 1+[(a:m-a)'; (m+a:2*m-a)'];

for k = 1:lyt
  B(k, pos(k)) = 1;
end
A = [sparse(A); B];
quadr = [wt, A];