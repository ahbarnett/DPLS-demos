function [A,An] = LapSLPmatrix(t,s,a) % single-layer kernel matrix & targ n-deriv
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
end

function i = diagind(A) % return indices of diagonal of square matrix
N = size(A,1); i = sub2ind(size(A), 1:N, 1:N);
end