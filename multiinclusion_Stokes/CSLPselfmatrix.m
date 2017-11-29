function S = CSLPselfmatrix(s,side) % complex SLP Kress-split Nystrom matrix
% s = src seg, even # nodes. side = 'i'/'e' int/ext case. Barnett 10/18/13.
% only correct for zero-mean densities.
N = numel(s.x);
d = repmat(s.x, [1 N]) - repmat(s.x.', [N 1]);  % C-# displacements mat, t-s
x = exp(1i*s.t); % unit circle nodes relative to which angles measured.
S = -log(d) + log(repmat(x,[1 N]) - repmat(x.',[N 1])); % NB circ angles
S(diagind(S)) = log(1i*x./s.xp);               % complex diagonal limit
% O(N^2) hack to remove 2pi phase jumps in S (assumes enough nodes for smooth):
for i=1:numel(S)-1, p = imag(S(i+1)-S(i)); S(i+1) = S(i+1) - 2i*pi*round(p/(2*pi)); end % (NB access matrix as 1d array!)
%figure; imagesc(real(S)); figure; imagesc(imag(S)); stop % check S mat smooth!

m = 1:N/2-1; if side=='e', Rjn = ifft([0 1./m 1/N 0*m]); % imag sign dep on side
else, Rjn = ifft([0 0*m 1/N 1./m(end:-1:1)]); end % cmplx Kress Rj(N/2)/4pi
%m = 1:N/2-1; Rjn = ifft([0 0*m 1/N 1./m(end:-1:1)]); Rjn = [Rjn(1) Rjn(end:-1:2)]; % flips order
S = S/N + circulant(Rjn); % incl 1/2pi SLP prefac. drops const part in imag
S = S .* repmat(s.sp.',[N 1]);  % include speed (2pi/N weights already in)
end