function [A,T] = DLPmatrix(t,s,mu,a) % double-layer 2D Stokes vel kernel matrix
% Returns 2N-by-2N matrix from src force vector to 2 flow components
% If detects self-int, does correct diagonal limit, no jump condition (PV int).
% t = target seg (x cols), s = src seg, a = optional transl of src seg.
% 2nd output is optional traction matrix, without self-eval. 2/1/14, 3/2/14
if nargin<4, a = 0; end
N = numel(s.x); M = numel(t.x);
% r = repmat(t.x, [1 N]) - repmat(s.x.' + a, [M 1]);    % C-# displacements mat
r = t.x*ones(1,N) - ones(M,1)*(s.x.'+a);
irr = 1./(conj(r).*r);    % 1/R^2
d1 = real(r); d2 = imag(r);
%rdotny = d1.*repmat(real(s.nx)', [M 1]) + d2.*repmat(imag(s.nx)', [M 1]);
rdotny = d1.*(ones(M,1)*real(s.nx)') + d2.*(ones(M,1)*imag(s.nx)');
rdotnir4 = rdotny.*(irr.*irr); if nargout<=1, clear rdotny; end
A12 = (1/pi)*d1.*d2.*rdotnir4;  % off diag vel block
A = [(1/pi)*d1.^2.*rdotnir4,   A12;                     % Ladyzhenzkaya
     A12,                      (1/pi)*d2.^2.*rdotnir4];
if numel(s.x)==numel(t.x) && max(abs(s.x+a-t.x))<1e-14
  c = -s.cur/2/pi;           % diagonal limit of Laplace DLP
  tx = 1i*s.nx; t1=real(tx); t2=imag(tx);     % tangent vectors on src curve
  A(sub2ind(size(A),1:N,1:N)) = c.*t1.^2;       % overwrite diags of 4 blocks:
  A(sub2ind(size(A),1+N:2*N,1:N)) = c.*t1.*t2;
  A(sub2ind(size(A),1:N,1+N:2*N)) = c.*t1.*t2;
  A(sub2ind(size(A),1+N:2*N,1+N:2*N)) = c.*t2.^2;
end
% A = A .* repmat([s.w(:)' s.w(:)'], [2*M 1]);            % quadr wei
A = A .* (ones(2*M,1)*[s.w(:)' s.w(:)']);
if nargout>1           % traction, my formula
  %nx1 = repmat(real(t.nx), [1 N]); nx2 = repmat(imag(t.nx), [1 N]);
   nx1 = real(t.nx)*ones(1,N); nx2 = imag(t.nx)*ones(1,N);
  rdotnx = d1.*nx1 + d2.*nx2;
%   ny1 = repmat(real(s.nx)', [M 1]); ny2 = repmat(imag(s.nx)', [M 1]);
  ny1 = ones(M,1)*real(s.nx)'; ny2 = ones(M,1)*imag(s.nx)';
  dx = rdotnx.*irr; dy = rdotny.*irr; dxdy = dx.*dy;
  R12 = d1.*d2.*irr; R = [d1.^2.*irr, R12; R12, d2.^2.*irr];
  nydotnx = nx1.*ny1 + nx2.*ny2;
  nyd=nydotnx.*irr - 8*dxdy;
  dxny1=dx.*ny1; dxny2=dx.*ny2; dyd1=dy.*d1; dyd2=dy.*d2;
  T = R.*[nyd nyd;nyd nyd] + [dxdy,0*dxdy;0*dxdy,dxdy];
    T = T + ([nx1.*ny1, nx1.*ny2; nx2.*ny1, nx2.*ny2] + ...
      [dxny1.*d1, dxny1.*d2; dxny2.*d1, dxny2.*d2] + ...
      [dyd1.*nx1, dyd1.*nx2; dyd2.*nx1, dyd2.*nx2]).*[irr irr;irr irr]; 
  
%   T = R.*kron(ones(2), nydotnx.*irr - 8*dxdy) + kron(eye(2), dxdy);
%   T = T + [nx1.*ny1.*irr, nx1.*ny2.*irr; nx2.*ny1.*irr, nx2.*ny2.*irr] + ...
%       kron(ones(2),dx.*irr) .* [ny1.*d1, ny1.*d2; ny2.*d1, ny2.*d2] + ...
%       kron(ones(2),dy.*irr) .* [d1.*nx1, d1.*nx2; d2.*nx1, d2.*nx2];
  T = (mu/pi) * T;                                        % prefac
  %T = T .* repmat([s.w(:)' s.w(:)'], [2*M 1]);            % quadr wei
  T = T.*(ones(2*M,1)*[s.w(:)' s.w(:)']); 
end
end