function [A,T] = SLPmatrix(t,s,mu,a) % single-layer 2D Stokes kernel vel matrix
% Returns 2N-by-2N matrix from src force vector to 2 flow component
% t = target seg (x cols), s = src seg, a = optional translation of src seg
% No option for self-int, gives Inf on diag.   3/2/14
% 2nd output is traction matrix, needs t.nx normal (C-#); no self-int either.
if nargin<4, a = 0; end
N = numel(s.x); M = numel(t.x);
% r = repmat(t.x, [1 N]) - repmat(s.x.' + a, [M 1]);    % C-# displacements mat
r = t.x*ones(1,N) - ones(M,1)*(s.x.' + a);
irr = 1./(conj(r).*r);    % 1/r^2
d1 = real(r); d2 = imag(r);
Ilogr = -log(abs(r));  % log(1/r) diag block
c = 1/(4*pi*mu);       % factor from Hsiao-Wendland book, Ladyzhenskaya
A12 = c*d1.*d2.*irr;   % off diag vel block
A = [c*(Ilogr + d1.^2.*irr), A12;                         % u_x
     A12,                c*(Ilogr + d2.^2.*irr)];         % u_y
% A = A .* repmat([s.w(:)' s.w(:)'], [2*M 1]);              % quadr wei
A = A .*(ones(2*M,1)*[s.w(:)' s.w(:)']);
if nargout>1           % traction (negative of DLP vel matrix w/ nx,ny swapped)
  %rdotn = d1.*repmat(real(t.nx), [1 N]) + d2.*repmat(imag(t.nx), [1 N]);
  rdotn = d1.*(real(t.nx)*ones(1,N)) + d2.*(imag(t.nx)*ones(1,N));
  rdotnir4 = rdotn.*(irr.*irr); clear rdotn
  A12 = -(1/pi)*d1.*d2.*rdotnir4;
  T = [-(1/pi)*d1.^2.*rdotnir4,   A12;                   % own derivation
     A12,                      -(1/pi)*d2.^2.*rdotnir4];
%   T = T .* repmat([s.w(:)' s.w(:)'], [2*M 1]);            % quadr wei
    T = T .*(ones(2*M,1)*[s.w(:)' s.w(:)']);     
end
end