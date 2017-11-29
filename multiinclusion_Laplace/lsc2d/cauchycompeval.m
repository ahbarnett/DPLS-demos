function [vc vcp] = cauchycompeval(x,s,vb,side,o)
% CAUCHYCOMPEVAL - global compensated evaluation int/exterior Cauchy integral
%
% v = cauchycompeval(x,s,vb,side) approximates the holomorphic function v
%  at targets x exterior to a closed curve s on which its complex boundary
%  values vb (v^+ for exterior case, v^- for interior) are given at nodes.
%  "side" specifies if the targets in interior or exterior to the curve.
%  For the exterior case, uniqueness of the holomorphic function is insured
%  by assuming the limit as |z|->infty is zero.
%  For both interior and exterior cases, this is done by evaluating the
%  Cauchy integral over the curve Gamma,
%
%      v(x) = (1/(2.pi.i)) integral_Gamma phi(y) / (x-y) dy,   where phi(y)=v(y)
%
%  using a barycentric-type formula accurate up to the boundary.
%
%  Warning!: the evaluation will not be any more accurate than the naive
%  trapezoid rule for generic complex data phi(y) in the above integral. High
%  accuracy is only achieved when phi lies in the appropriate Hardy space, that
%  is, when phi is the boundary data of some function holomorphic either inside,
%  or holomorphic outside with correct decay condition.
%
% [v vp] = cauchycompeval(x,s,vb,side) also returns v's complex derivative
%
% [v ...] = cauchycompeval(x,s,vb,side,opts) controls various options.
%
% For both values and derivatives, true barycentric forms are used that ensure
% machine accuracy for arbitrarily small (or zero) target-node distances.
%
% Inputs:
% x = row or col vec of M targets, as points in complex plane
% s = closed curve quadrature struct containing N nodes s.x, and s.w
%     "speed weights", s.nx unit normals, s.a interior pt far from the bdry.
%     All vectors are represented as complex numbers.
% vb = row or col vec of N boundary values of holomorphic function v
% side = 'i' or 'e' appropriate for if targets interior or exterior to curve.
% opts = options struct controlling the following:
%        opts.delta : override default setting of delta (min dist O(N) derivs)
%                     (set delta=0 to use naive non-bary derivative formula)
% Outputs:
% v  = row vec approximating the homolorphic function v at the M targets
% vp = (optional) row vec of complex first derivative v' at the M targets
%
% Notes; 1) algorithm is that of Ioakimidis et al BIT 1991 for interior, and,
% for the exterior, a modified version using 1/(z-a) in place of the function 1.
% (The exterior case could be done slightly more simply using Helsing's
% proposal from (26)-(27) of his 2008 JCP paper).
% For the derivative, the formula is mathematically the derivative of the
% barycentric formula of Schneider-Werner 1986 described in Berrut et al 2005,
% but using the complex quadrature weights instead of the barycentric weights.
% However, since this derivative is not itself a true barycentric form, a hack
% is needed to compute the difference v_j-v(x) in a form where roundoff error
% cancels correctly. The exterior case is slightly more intricate. When
% a target coincides with a node j, the true value v_j (or Schneider-Werner
% formula for v'(x_j)) is used.
% 2) In order to vectorize in both the node and target indices, we use a lot
% of RAM, O(NM), which is not really necessary. For now, the user should group
% their targets into smaller blocks done with separate calls if the RAM use
% is too high.
%
% Alex Barnett 10/22/13 based on cauchycompevalint, pole code in lapDevalclose.m
% 10/23/13 node-targ coincidences fixed, exterior true barycentric discovered.
%
% To do: code up and compare to Helsing's exterior formula.
if nargin<5, o = []; end

if ~isfield(o,'delta'), o.delta = 1e-2*0; end
% dist param where $O(N^2.M) deriv bary form switched on.
% Roughly deriv errors are then limited to emach/mindist
% The only reason to decrease mindist is if lots of v close
% nodes on the curve with lots of close target points.

cw = s.cw;                    % complex speed weights
N = numel(s.x); M = numel(x);

if nargout==1  % no deriv wanted... (note sum along 1-axis faster than 2-axis)
  comp = repmat(cw(:), [1 M]) ./ (repmat(s.x(:),[1 M]) - repmat(x(:).',[N 1]));
  if side=='e', pcomp = comp .* repmat(1./(s.x(:)-s.a), [1 M]);
  else pcomp = comp; end  % pcomp are weights and bary poles appearing in J0
  I0 = sum(repmat(vb(:),[1 M]).*comp); J0 = sum(pcomp); % Ioakimidis notation
  vc = I0./J0;                        % bary form
  if side=='e', vc = vc./(x(:).'-s.a); end         % correct w/ pole
  [jj ii] = ind2sub(size(comp),find(~isfinite(comp))); % node-targ coincidences
  for l=1:numel(jj), vc(ii(l)) = vb(jj(l)); end % replace each hit w/ corresp vb
  
else           % 1st deriv wanted...
  invd = 1./(repmat(s.x(:),[1 M]) - repmat(x(:).',[N 1])); % 1/displacement mat
  comp = repmat(cw(:), [1 M]) .* invd;
  if side=='e', pcomp = comp .* repmat(1./(s.x(:)-s.a), [1 M]);
  else pcomp = comp; end  % pcomp are weights and bary poles appearing in J0
  I0 = sum(repmat(vb(:),[1 M]).*comp); J0 = sum(pcomp);
  if side=='e', prefac = 1./(x(:).'- s.a); else prefac = 1; end
  vc = prefac .* I0./J0;                        % bary form (poss overall pole)
  dv = repmat(vb(:),[1 M]) - repmat(vc(:).',[N 1]); % v value diff mat
  [jj ii] = ind2sub(size(invd),find(abs(invd) > 1/o.delta)); % bad pairs indices
  if side=='e' % exterior:
    for l=1:numel(jj), j=jj(l); i=ii(l);p = sum(comp(:,i).*(vb(j)./(s.x(:)-s.a)-vb(:)./(s.x(j)-s.a))) / sum(pcomp(:,i));dv(j,i) = prefac(i) * ((s.x(j)-s.a)*p + (x(i)-s.x(j))*vb(j));end % pole-corrected for dv, gives bary stability for close node-targ pairs
  else  % interior:
    for l=1:numel(jj), j=jj(l); i=ii(l); % loop over node-targ pairs too close
      dv(j,i) = sum(comp(:,i).*(vb(j)-vb(:))) / sum(pcomp(:,i));
    end % bary for dv, gives bary stability for close node-targ pairs
  end
  vcp = prefac .* sum(dv.*comp.*invd) ./ J0; % bary form for deriv
  [jj ii] = ind2sub(size(comp),find(~isfinite(comp))); % node-targ coincidences
  for l=1:numel(jj), j=jj(l); i=ii(l); % loop over hitting node-targ pairs
    vc(i) = vb(j);                     % replace each hit w/ corresp vb
    notj = [1:j-1, j+1:N];             % Schneider-Werner form for deriv @ node:
    vcp(i) = -sum(cw(notj).*(vb(j)-vb(notj))./(s.x(j)-s.x(notj)))/cw(j);
  end
end
