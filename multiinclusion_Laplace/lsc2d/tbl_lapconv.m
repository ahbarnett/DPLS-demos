function tbl_lapconv
% edit of testlapSDevalclose for lsc2d paper table. Barnett 11/11/13
% Now uses external quadr(). 10/8/14

a = .3; w = 5;         % smooth wobbly radial shape params
R = @(t) 1 + a*cos(w*t); Rp = @(t) -w*a*sin(w*t); Rpp = @(t) -w*w*a*cos(w*t);
s.Z = @(t) R(t).*exp(1i*t); s.Zp = @(t) (Rp(t) + 1i*R(t)).*exp(1i*t);
s.Zpp = @(t) (Rpp(t) + 2i*Rp(t) - R(t)).*exp(1i*t);
s.inside = @(z) abs(z)<R(angle(z)); s.a = -0.1;  % interior pt far from bdry

for expt = {'di','de','si','se'}  % =========== four cases
  lp = expt{1}(1), side = expt{1}(2) % choose expt: SLP/DLP and int/ext (4 cases)
if side=='i'   % interior
  p = []; p.x = .2+.3i; po=p; po.x = 0;  % test pt & offset-computing pt for Neu
  fholom = @(z) exp(1i*(z+1));        % holomorphic exact interior soln
  fpholom = @(z) 1i*exp(1i*(z+1));        % its complex derivative
  %a = [1.5+1.5i, -0.25+1.5i, -.5-1.5i];   % ...or, Helsing 2008 data Sec 10.1
  %fholom = @(z) 1./(z-a(1))+1./(z-a(2))+1./(z-a(3));  % 3 poles
  %fpholom = @(z) -1./(z-a(1)).^2-1./(z-a(2)).^2-1./(z-a(3)).^2; % (Hels done)
else       % exterior (note has O(1/|z|) decay so not general Dir BVP soln)
  p = []; p.x = 1.2+1i;
  a=.1+.3i; fholom = @(z) 1./(z-a);   % holomorphic exact exterior soln
  fpholom = @(z) -1./(z-a).^2;        % its complex derivative
end
f = @(z) real(fholom(z)); imf = @(z) imag(fholom(z)); % real part (true f), etc
fx = @(z) real(fpholom(z)); fy = @(z) -imag(fpholom(z)); % note sign!

nx = 300; gx = ((1:nx)/nx*2-1)*1.5; gy = gx; % plotting grid fror BVP test
[xx yy] = meshgrid(gx,gy); zz = (xx+1i*yy); clear xx yy
fz=f(zz); fxz=fx(zz); fyz=fy(zz); % store exact fields on grid

ii = s.inside(zz); if side=='e', ii = ~ii; end % ii = indices to eval at

Ns = 50:50:250;
es = nan(2,numel(Ns));
for i=1:numel(Ns), N = Ns(i); % ------------------- convergence loop
  s = quadr(s,N); % set up nodes on bdry (note outwards normal)

% BIE & solve: use jump relation to build A which maps density to values on s...
uoff = 0;                            % no default offset
if lp=='d' % DLP test via Dirichlet BVP
  rhs = [f(s.x)];                    % col vec of Dirichlet data
  if side=='i'
    A = -eye(N)/2 + DLPmatrix(s,s);  % full rank
  else
    A = eye(N)/2 + DLPmatrix(s,s);   % has rank-1 nullspace, ok for dense solve
  end
  tau = A \ rhs;                     % solve coefficients
  u = DLPmatrix(p,s) * tau;          % check eval @ test pt
elseif lp=='s'  % SLP test via Neu BVP
  rhs = [real(s.nx).*fx(s.x)+imag(s.nx).*fy(s.x)];  % col vec of Neu data 
  if side=='i'  % nonunique so will need to offset by const using the origin...
    [~,A] = SLPmatrix(s,s); A = A + eye(N)/2; % nullity 1 but consistent
  else
    [~,A] = SLPmatrix(s,s); A = A - eye(N)/2; % full rank
  end
  tau = A \ rhs;                     % solve coefficients
  u = SLPmatrix(p,s) * tau;          % check eval @ test pt
  if side=='i', uoff = f(po.x) - SLPmatrix(po,s) * tau; % measure offset
  %else, fprintf('tau integral = %.3g (zero for ext Neu BVP)\n',sum(tau.*s.w))
  end
end
%fprintf('native u err @ dist pt (ie BVP soln err) = %.3g\n', abs(u+uoff-f(p.x)))

z = s.x(ceil(N/3)) + 1e-12;   % check close eval at one point near to or on a node...
if lp=='s', [u ux uy] = lapSevalclose(z,s,tau,side);
  else, [u ux uy] = lapDevalclose(z,s,tau,side); end
%fprintf('u,ux,uy near pt err = %.3g,%.3g,%.3g\n',u+uoff-f(z),ux-fx(z),uy-fy(z));

%fprintf('close eval of pots and flds... ');
  uc = nan*zz; uxc = uc; uyc = uc; tic % arrays for answers
  if lp=='d',[uc(ii) uxc(ii) uyc(ii) do] = lapDevalclose(zz(ii(:)),s,tau,side);
  else, [uc(ii) uxc(ii) uyc(ii) do] = lapSevalclose(zz(ii(:)),s,tau,side); end
  t=toc; %fprintf('done in %.3g s (%.3g src-targ pairs/s)\n',t, N*numel(p.x)/t)
  
% report error sup norm outputs and maybe plot...
es(1,i) = max(abs(uc(:)+uoff-fz(:)));
%max(abs(uc(:)))  % sup norm of solution
es(2,i) = max(abs(uxc(:)-fxz(:)+1i*uyc(:)-1i*fyz(:)));
fprintf('N=%d:     \t%.2g    \t%.2g\n', N, es(1,i), es(2,i))
end % -------------------------------- end conv loop
end % ======= expt loop
%keyboard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function h = showseg(s, uc) % plot a segment & possibly its nei copies
if nargin<2, uc.nei = 0; uc.d = 0; end % dummy uc
if uc.nei>0, uc.nei = 1; end % don't plot a ridiculous # of copies
for i=-uc.nei:uc.nei
  plot(s.x+uc.d*i, 'b.-'); axis equal xy;
  l=0.3; hold on; plot([s.x+uc.d*i, s.x+l*s.nx+uc.d*i].', 'k-');
end

function [A An] = DLPmatrix(t,s,a) % double-layer kernel matrix & targ n-deriv
% t = target seg (x,nx cols), s = src seg, a = optional translation of src seg
% No jump included on self-interaction (ie principal-value integral)
if nargin<3, a = 0; end
N = numel(s.x); M = numel(t.x);
d = repmat(t.x, [1 N]) - repmat(s.x.' + a, [M 1]);    % C-# displacements mat
ny = repmat(s.nx.', [M 1]);      % identical rows given by src normals
A = (1/2/pi) * real(ny./d);      % complex form of dipole
if numel(s.x)==numel(t.x) & max(abs(s.x+a-t.x))<1e-14
  A(diagind(A)) = -s.cur/4/pi; end           % self? diagonal term for Laplace
A = A .* repmat(s.w(:)', [M 1]);
if nargout>1 % deriv of double-layer. Not correct for self-interaction.
  csry = conj(ny).*d;              % (cos phi + i sin phi).r
  nx = repmat(t.nx, [1 N]);        % identical cols given by target normals
  csrx = conj(nx).*d;              % (cos th + i sin th).r
  r = abs(d);    % dist matrix R^{MxN}
  An = -real(csry.*csrx)./(r.^2.^2)/(2*pi);
  An = An .* repmat(s.w(:)', [M 1]);
end

function [A An] = SLPmatrix(t,s,a) % single-layer kernel matrix & targ n-deriv
% t = target seg (x,nx cols), s = src seg, a = optional translation of src seg
% No jump included on self-interaction of derivative (ie PV integral).
if nargin<3, a = 0; end
N = numel(s.x); M = numel(t.x);
d = repmat(t.x, [1 N]) - repmat(s.x.' + a, [M 1]);    % C-# displacements mat
A = -(1/2/pi) * log(abs(d)) .* repmat(s.w(:)', [M 1]);  % infty for self-int
if nargout==2                      % apply D^T
  nx = repmat(-t.nx, [1 N]);       % identical cols given by -targ normals
  An = (1/2/pi) * real(nx./d);     % complex form of dipole. Really A1 is An
  if numel(s.x)==numel(t.x) & max(abs(s.x+a-t.x))<1e-14
    An(diagind(An)) = -s.cur/4/pi; end  % self? diagonal term for Laplace
  An = An .* repmat(s.w(:)', [M 1]);
end

function i = diagind(A) % return indices of diagonal of square matrix
N = size(A,1); i = sub2ind(size(A), 1:N, 1:N);
