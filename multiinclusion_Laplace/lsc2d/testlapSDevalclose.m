function testlapSDevalclose
% Solve 2D Dirichlet Laplace BVPs to test global close-evaluation formulae
% of Barnett which generalize Helsing/Ioakimidis.
% Currently DLP and SLP (lp='s','d'), int and ext cases (side='i','e').
% We use Ch.6 of LIE Kress book for formulation and uniqueness of the 4 BVPs,
% except we only consider u_infty=0 for ext Dir BVP (and don't modify DLP here).
% For ext Neu BVP only the zero-mean SLP is tested for close-eval.
% Needs: lapDevalclose.m, lapSevalclose.m, quadr.m, and whatever they needs.
% Alex Barnett for Shravan + Wu 10/18/13-11/4/13, based on 4/20/13.

clear; v = 1; der = 1; % verbosity (v=0,1,2); der=0 (no derivs), 1 (test derivs)
lp = 's'; side = 'e';  % choose expt: SLP/DLP and int/ext (one of 4 cases)
a = .3; w = 5;         % smooth wobbly radial shape params
R = @(t) 1 + a*cos(w*t); Rp = @(t) -w*a*sin(w*t); Rpp = @(t) -w*w*a*cos(w*t);
s.Z = @(t) R(t).*exp(1i*t); s.Zp = @(t) (Rp(t) + 1i*R(t)).*exp(1i*t);
s.Zpp = @(t) (Rpp(t) + 2i*Rp(t) - R(t)).*exp(1i*t);
s.inside = @(z) abs(z)<R(angle(z)); s.a = -.1+.2i;  % interior pt far from bdry
N = 220; s = quadr(s,N); % set up nodes on bdry (note outwards normal)

if side=='i'   % interior
  p = []; p.x = .2+.3i; po=p; po.x = 0;  % test pt & offset-computing pt for Neu
  fholom = @(z) exp(1i*(z+1));        % holomorphic exact interior soln
  fpholom = @(z) 1i*exp(1i*(z+1));        % its complex derivative
else       % exterior (note has O(1/|z|) decay so not general Dir BVP soln)
  p = []; p.x = 1.2+1i;
  a=.1+.3i; fholom = @(z) 1./(z-a);   % holomorphic exact exterior soln
  fpholom = @(z) -1./(z-a).^2;        % its complex derivative
end
f = @(z) real(fholom(z)); imf = @(z) imag(fholom(z)); % real part (true f), etc
fx = @(z) real(fpholom(z)); fy = @(z) -imag(fpholom(z)); % note sign!

% Initial test given density, if value and deriv agrees in the far-field where
% standard trapezoid rule is accurate. Can test if non-zero mean SLP is good:
tau = exp(sin(4*s.t+1));  % generic smooth density, non-zero mean
p.nx = exp(1i*pi/3); % normal direction, at target pt p chosen above
if lp=='d', [A An] = DLPmatrix(p,s); [u ux uy] = lapDevalclose(p.x,s,tau,side);
else, [A An] = SLPmatrix(p,s); [u ux uy] = lapSevalclose(p.x,s,tau,side); end
u0 = A*tau; un0 = An*tau;
un = real(p.nx)*ux + imag(p.nx)*uy; % close-eval target normal deriv
fprintf('close-eval minus native quadr @ dist pt: val = %.3g, nder = %.3g\n', u-u0, un-un0)

nx = 300; gx = ((1:nx)/nx*2-1)*1.5; gy = gx; % plotting grid fror BVP test
[xx yy] = meshgrid(gx,gy); zz = (xx+1i*yy); clear xx yy
fz=f(zz); fxz=fx(zz); fyz=fy(zz); % store exact fields on grid
if v>1, figure; imagesc(gx,gy,fz); colormap(jet(256)); colorbar;
  caxis([-2 2]); hold on; showseg(s); title('exact solution = f'); end

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
  else fprintf('tau integral = %.3g (zero for ext Neu BVP)\n',sum(tau.*s.w)),end
end
fprintf('native u err @ dist pt (ie BVP soln err) = %.3g\n', abs(u+uoff-f(p.x)))

ii = s.inside(zz); if side=='e', ii = ~ii; end  % ii = indices to eval at
p.x = zz(ii(:));  % eval pts only on one side
ug = nan*zz;
if lp=='d', ug(ii) = DLPmatrix(p,s) * tau;  % native eval u on grid (uses RAM)
else, ug(ii) = SLPmatrix(p,s) * tau; end
ug = ug + uoff; ug = reshape(ug, size(zz)); % include non-uniqueness offset
if v, figure; imagesc(gx,gy,log10(abs(ug-f(zz)))); colormap(jet(256));
  colorbar; caxis([-16 0]); title('log_{10} error in u via native bdry quadr');
  hold on; showseg(s); end
  
z = s.x(165) + 1e-12;   % check close eval at one point near to or on a node...
if der, if lp=='s', [u ux uy] = lapSevalclose(z,s,tau,side);
  else, [u ux uy] = lapDevalclose(z,s,tau,side); end
fprintf('u,ux,uy near pt err = %.3g,%.3g,%.3g\n',u+uoff-f(z),ux-fx(z),uy-fy(z));
else, if lp=='s', u = lapSevalclose(z,s,tau,side);
  else, u = lapDevalclose(z,s,tau,side); end
  fprintf('u near pt err = %.3g\n',u+uoff-f(z)); end

if der==0, fprintf('close eval of pots... ');  % close eval on the grid
  uc=nan*zz; tic, if lp=='d',[uc(ii)] = lapDevalclose(zz(ii(:)), s, tau, side);
  else, [uc(ii)] = lapSevalclose(zz(ii(:)), s, tau, side); end
  t=toc; fprintf('done in %.3g s (%.3g src-targ pairs/s)\n',t, N*numel(p.x)/t)
else, fprintf('close eval of pots and flds... ');
  uc = nan*zz; uxc = uc; uyc = uc; tic % arrays for answers
  if lp=='d',[uc(ii) uxc(ii) uyc(ii) do] = lapDevalclose(zz(ii(:)),s,tau,side);
  else, [uc(ii) uxc(ii) uyc(ii) do] = lapSevalclose(zz(ii(:)),s,tau,side); end
  t=toc; fprintf('done in %.3g s (%.3g src-targ pairs/s)\n',t, N*numel(p.x)/t)
end
  
if v>1 & der % test Helsing step 1 reprod holomorphic f: off by imag const
  figure; subplot(2,1,1); plot(real(do.vb)+uoff-f(s.x), '-'); title('Re [vb-f]'); subplot(2,1,2); plot(imag(do.vb)-imf(s.x), '-'); title('Im [vb-f]'); end

% report error sup norm outputs and maybe plot...
fprintf('u:   max = %.3g, max close-eval err = %.3g\n', max(abs(fz(ii(:)))), max(abs(uc(:)+uoff-fz(:))))
if der
  fprintf('u_x: max = %.3g, max close-eval err = %.3g\n', max(abs(fxz(ii(:)))), max(abs(uxc(:)-fxz(:))))
  fprintf('u_y: max = %.3g, max close-eval err = %.3g\n', max(abs(fyz(ii(:)))), max(abs(uyc(:)-fyz(:))))
end
if v, figure; imagesc(gx,gy,log10(abs(uc+uoff-fz))); colormap(jet(256));
  colorbar; title('log_{10} error in u'); hold on; showseg(s);
  if der, figure; imagesc(gx,gy,log10(abs(uxc-fx(zz)))); colormap(jet(256));
    colorbar; title('log_{10} error in u_x'); hold on; showseg(s);
    figure; imagesc(gx,gy,log10(abs(uyc-fy(zz)))); colormap(jet(256));
    colorbar; title('log_{10} error in u_y'); hold on; showseg(s); end
end
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
