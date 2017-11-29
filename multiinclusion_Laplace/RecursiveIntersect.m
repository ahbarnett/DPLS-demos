function t=RecursiveIntersect(t,n)

t0M=floor(length(t)/2/n);

t0=t(1:n*t0M,:);
t1=t(n*t0M+1:end,:);

t0=IslandIntersect(t0,t1,n);
t1=IslandIntersect(t1,t0,n);


if length(t0)/n>1
    t0=RecursiveIntersect(t0,n);
end
    
if length(t1)/n>1
   t1=RecursiveIntersect(t1,n);
end
t=[t0;t1];