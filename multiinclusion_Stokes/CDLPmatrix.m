function CDLP=CDLPmatrix(s)

CDLP=cell(1,2);

N=numel(s.x);
Sj=ones(N,1)*s.x.';
Wj=ones(N,1)*s.cw.';
CDLP{2}=Wj./(Sj-Sj.');
CDLP{2}(diag(true(N,1)))=0;

beta = 2.2;
N = ceil(beta*numel(s.x)/2)*2;
sf.x = fftinterp(s.x,N); sf = quadr(sf);
Sj=ones(N,1)*sf.x.';
Wj=ones(N,1)*sf.cw.';
CDLP{1}=Wj./(Sj-Sj.');
CDLP{1}(diag(true(N,1)))=0;
end

function g = fftinterp(f,N)
% FFTINTERP - resample periodically sampled function onto finer grid
%
% g = fftinterp(f,N)
% inputs:  f - (row or column) vector length n of samples
%          N - desired output number of samples, must be >= n
% outputs: g - vector length N of interpolant, (row or col as f was)
% Note on phasing: the output and input grid first entry align.
% Barnett 9/5/14.  To do: downsample case N<n
n = numel(f);
if N==n, g = f; return; end
if mod(N,2)~=0 || mod(n,2)~=0, warning('N and n must be even'); end
F = fft(f(:).');    % row vector
g = ifft([F(1:n/2) F(n/2+1)/2 zeros(1,N-n-1) F(n/2+1)/2 F(n/2+2:end)]);
g = g*(N/n);   % factor from the ifft
if size(f,1)>size(f,2), g = g(:); end % make col vector
end