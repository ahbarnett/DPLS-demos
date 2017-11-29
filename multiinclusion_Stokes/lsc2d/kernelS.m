function G = kernelS(X)
%KERNELS returns the integration kernel for the single layer potential
%over the contour X. X should be a single column.
%

persistent qw;

if(isempty(qw))
%    qw = quadratureS(length(X)/4, 8, 2*pi);  %Alpert's 8th order quadrature for singular integrals
   qw = quadratureS(length(X)/4, 16, 2*pi);  %Alpert's 16th order quadrature for singular integrals
end

N = length(X)/2; X = [X(1:N), X(N+1:2*N)];

w = qw(:,1);  A = qw(:, 2:N+1);
G = zeros(2*N,2*N);  
for j =1:N
    ind = 1+ mod(j-1 + (0:N-1), N); Xin = [A*X(ind,1), A*X(ind,2)]; 
    
    rho = (sum((Xin - repmat(X(j,:), length(Xin), 1)).^2, 2));
    br = 1./rho;

    x = X(j,:); 
    y = Xin;    
    LogPart = -0.5* (w.*log(rho) )'*A;  
    G(j, ind) = (LogPart + (w.*( (x(1) -y(:, 1)).^2.*br))'*A);  
    G(j, N+ind) = ((w.*(x(1) -y(:, 1)).*(x(2) - y(:, 2)).*br)'*A);
    G(j+N, 1:N) = G(j, N+1:2*N);
    G(j+N, N+ind) = (LogPart + (w.*( (x(2) -y(:, 2)).^2.*br))'*A);
end

sa = sqrt(DmFT(X(:,1),1).^2 + DmFT(X(:,2),1).^2);
G = G.*[repmat(sa', 2*N, 1), repmat(sa', 2*N, 1)];
clear qw

%   for rr = 0
%     
%   n = 64; nv = 1;
%   g = (1:n)'*2*pi/n; R = 10^rr;
%   X = R*[cos(g) sin(g)];
%   
%   u = zeros(2*n,1);
%   sig = ones(n,1)*R;
%   
%   prams.kappa = R^3; 
%   options.useFMM = 0;
%   
%   viscCont = 1; 
% % On a circle with this kappa and sig , f = -n.
%   
%   [F,FX] = InteractionForce(X(:),u,sig,prams,options,viscCont);
% 
%   f = -X/R;
%   Fex = X/4; % r \otimes r part
%   lambda = abs(8*[1 1:n/2-1 -n/2:-1]');
%   lambda = [lambda lambda];
%   inv_lambda = 1./lambda; inv_lambda(1,:) = 0; 
%   Fex = ifft(inv_lambda.*fft(f)); %log part
%    
%   e1 = Fex-F; e2 = sqrt(dot(Fex,Fex,2)./dot(F,F,2));
%   disp([rr max(dot(e1,e1,2)) mean(e2)]);
%   disp('----------------------------------');
%   end
%   %quiver(X(:,1),X(:,2),e2(:,1),e2(:,2));
