function df = DmFT(f, m, kIn, saIn)
%DmFT(f,m) Computes the mth derivative of a periodic function f using
%fourier transform. By calling DmFT(f,m,K,Sa) one can overrride the
%differentiation coefficients with the desired ones by passing extra
%arguments K and the jacobian Sa.
%
%EXAMPLE:
%  n = 64;
%  X = boundary(n,'nv',1); X = reshape(X,n,2);
%  DX  = DmFT(X,1);
%  DDX = DmFT(X,2);
%
%  sa = sqrt(dot(DX,DX,2));
%  tang = DX./[sa sa];
%  kap = (DX(:,1).*DDX(:,2) - DX(:,2).*DDX(:,1))./sa.^3;
% 
%  viewer(X(:),tang(:)); 
%

  persistent k sa


  if(nargin == 0), testDmFT(); return; end
  if(nargin == 1 || isempty(m)), m=1 ;end
  if(nargin >= 3), k = kIn; end
  if(nargin == 4), sa = saIn; end

  if(isempty(k))
    N = size(f,1);
    if(mod(N,2)==1),warning('DmFT:size',['Are you sure about using FFT ' ...
                          'for differentiation of a vector \n with odd ' ...
                          'length?'])
    end
    k = 1i*[(0:N/2-1)  0 (-N/2+1:-1)]';
  end

  if(isempty(sa)), sa = ones(size(f,1),1);end

  df = f;

  if(~isempty(f))
    [n, col] = size(f);
    
    if(size(sa,2)~=col)
      S = repmat(sa,1,col);
    else
      S = sa;  
    end
    
    K = repmat(k,1,col);
    
    for j = 1:m
      df = S.*real(ifft(K.*fft(df)));
    end
  end
  clear k sa

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function testDmFT()
% Source: Spectral Methods in Matlab, Nick Trefethen.

  disp('Running DmFT.m test:')
% Differentiation of a hat function:
  N = 24; h = 2*pi/N; x = h*(1:N)';
  v1 = max(0,1-abs(x-pi)/2);
  w = DmFT(v1,1,[],[]); clf
  subplot(3,2,1), plot(x,v1,'.-','markersize',13)
  axis([0 2*pi -.5 1.5]), grid on, title('function')
  subplot(3,2,2), plot(x,w,'.-','markersize',13)
  axis([0 2*pi -1 1]), grid on, title('spectral derivative')
  
% Differentiation of exp(sin(x)):
  v2 = exp(sin(x)); vprime = cos(x).*v2;
  w = DmFT(v2,1,[],[]);
  error = norm(w-vprime,inf);
  subplot(3,2,3), plot(x,v2,'.-','markersize',13)
  axis([0 2*pi 0 3]), grid on
  subplot(3,2,4), plot(x,w,'.-','markersize',13)
  axis([0 2*pi -2 2]), grid on
  text(2.2,1.4,['max error = ' num2str(error)])
  
  disp('-----------------------');
  disp('exp(sin(x)) derivative:');
  disp('-----------------------');
  for N=[8 16 32 64 128]
    h = 2*pi/N; x = h*(1:N)';
    v = exp(sin(x)); vprime = cos(x).*v;
    w = DmFT(v,1,[],[]);
    error = norm(w-vprime,inf);
    fprintf('  %3.0f points: |e|_inf = %2.4e\n',N, error);
  end
  
  disp('-------------------------------------------------------------------------');
  disp('cos(x)sin(x)+sin(x)cos(10x)+sin(20x)cos(13x) first and second derivative:');
  disp('-------------------------------------------------------------------------');
  for N=10:10:120
    h = 2*pi/N; x = h*(1:N)';
    v = cos(x).*sin(x)+sin(x).*cos(10*x)+sin(20*x).*cos(13*x);
    vprime = cos(2*x)-9/2*cos(9*x)+11/2*cos(11*x)+7/2*cos(7*x)+33/2*cos(33*x);
    v2prime = -4*cos(x).*sin(x) - 101*sin(x).*cos(10*x) - 20*cos(x).* ...
              sin(10*x)-569*sin(20*x).* cos(13*x)-520*cos(20*x).*sin(13*x);

    w = DmFT(v,1,[],[]);
    w2 = DmFT(v,2,[],[]);
    error = norm(w-vprime,inf);
    error2 = norm(w2-v2prime,inf);
    fprintf('  %3.0f points: |e|_inf = %2.4e,\t |e''|_inf = %2.4e\n',N, ...
            error,error2);
  end
  
  disp('-------------------------------------------------------------------------');
  fprintf(' m = 0 : %2.4e\n', norm(v-DmFT(v,0),inf));
  fprintf(' k = 0 : %2.4e\n', norm(DmFT(v,1,zeros(size(v))),inf));
  fprintf(' k = 1 : %2.4e\n', norm(v-DmFT(v,1,ones(size(v))),inf));
  fprintf('sa = 0 : %2.4e\n', norm(DmFT(v,1,[],zeros(size(v))),inf));
  fprintf('sa = 2 : %2.4e\n', norm(2*vprime-DmFT(v,1,[],2*ones(size(v))),inf));
  disp('-------------------------------------------------------------------------');