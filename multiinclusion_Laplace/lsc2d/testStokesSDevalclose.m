function testStokesSDevalclose(N,side,lp)
% testStokesSDevalclose(N,side,lp) produces a figure and error text data,
% testing a Stokes layer potential evaluation scheme.
%
% Input: N = # quadrature points on curve
%        side = 'i' for interior, 'e' for exterior
%        lp = 's' for single, 'd' for double
%             (or for BVP type:  's' for Neumann, 'd' for Dirichlet)
%
% Also see: StokesScloseeval, StokesDcloseeval, tbl_Stokesconv.m
%
% (c) Bowei Wu, Sept 2014. Barnett tweaked interface of Stokes{S,D}closeeval.

% close all
v = 1;  % 0=suppress figure outputs, 1=output figures, 2=plot only density functions
if nargin<3, lp = 's'; end     % note default is SLP
if nargin<2, side = 'e'; end
if nargin<1, N = 250; end

%% 1 generate the exact solution

% 1.1 set up domain: starfish domain
a = .3; w = 5;         % smooth wobbly radial shape params
R = @(t) 1 + a*cos(w*t); Rp = @(t) -w*a*sin(w*t); Rpp = @(t) -w*w*a*cos(w*t);
s.Z = @(t) R(t).*exp(1i*t); s.Zp = @(t) (Rp(t) + 1i*R(t)).*exp(1i*t);
s.Zpp = @(t) (Rpp(t) + 2i*Rp(t) - R(t)).*exp(1i*t);
s.inside = @(z) abs(z)<R(angle(z));
if side == 'e'
    s.inside = @(z) abs(z)<=R(angle(z));
end
s.a = -0.1;  % interior pt far from bdry
s = quadr(s, N);

% 1.2 set up locations of point forces
alpha = linspace(-pi/4,pi/4,5);
y_force = [];                             
if side == 'e'
    y_force.x = .3*exp(1i*alpha).';   % location of the point forces (exterior)
else
    y_force.x = 2*exp(1i*alpha).'; % location of the point forces (interior)
end
% pt_force = [[1;1;0;1;0];[0;0;1;0;1]];  % magnitude & direction of the point forces
pt_force = [[1;1;0;1;0];[1;0;-1;0;1]];  % magnitude & direction of the point forces (sample 2)

% 1.3.1 set up solution meshgrid
nx = 150; gx = ((1:nx)/nx*2-1)*1.5; gy = gx; % plotting grid
[xx yy] = meshgrid(gx,gy); zz = (xx+1i*yy); clear xx yy

% 1.3.2 use stokeslet to eval exact soln
fhom = nan*zz; % exact soln
p = [];
ii = s.inside(zz); if side=='e', ii = ~ii; end  % ii = indices to eval at
p.x = zz(ii(:));  % eval pts only on one side
A = stokesletmatrix(p,y_force);
f_temp = A*pt_force;
fhom(ii) = f_temp(1:end/2) + 1i*f_temp(end/2+1:end); % the exact soln (complex form)

% plot the exact soluntion and the point forces
if v==1
    figure; streamslice(gx,gy,real(fhom),imag(fhom));
    hold on; t = (1:N)'/N*2*pi; plot(real(s.Z(t)),imag(s.Z(t)),'r');
%     title('Exact Soln')
    plot(y_force.x,'.','MarkerSize',10,'LineWidth',10)
    quiver(real(y_force.x),imag(y_force.x),.2*pt_force(1:end/2),.2*pt_force(end/2+1:end))
    axis equal tight
end

%% 2 solve the BVP

% 2.1 get bdry condition
f = stokesletmatrix(s,y_force)*pt_force; % Dirichlet bdry condition
fp = stokesletpmatrix(s,y_force)*pt_force; % Neumann bdry condition

%'sup norm of bdry data f, fp:' max(abs(f(:))), max(abs(fp(:)))

% 2.2 solve for tau
warning('off','MATLAB:nearlySingularMatrix')
if lp == 's'
    A = SLPmatrixp(s,s); % normal deriv of SLP
    if side == 'e'
        tau = (-eye(size(A))/2 + A) \ fp;      % operator has rank-1 nullspace
    elseif side == 'i'
        tau = (eye(size(A))/2 + A) \ fp;       % operator has rank-3 nullspace
        % but is consistent, could easily fix
    end
elseif lp == 'd'
    if side == 'e'
        A = DLPmatrix(s,s)+kernelS([real(s.x);imag(s.x)]);   % operator has full rank, 
  % kernelS is SLP matrix filled by Alpert's singular quad
%         A = DLPmatrix(s,s);     % operator has rank-3 nullspace
        tau = (eye(size(A))/2 + A)\f; 
    elseif side == 'i'
        A = DLPmatrix(s,s); %+kernelS([real(s.x);imag(s.x)]); % for testing alpert errs    % operator has rank-1 nullspace
        tau = (-eye(size(A))/2 + A)\f;
    end
end
%     if v, figure, plot(tau), title('density'), end % plot density

% 2.3 eval SLP/DLP on the same meshgrid
uc = nan*zz; % layer potential
tau = tau(1:end/2)+1i*tau(end/2+1:end);
if lp == 's'
    uc(ii) = StokesScloseeval(p.x, s, tau, side);
    if side=='i'   % correct the rigid body (rank 3) nullspace problem -Alex
      z0 = [0;1];  % two test pts inside
      u0 = StokesScloseeval(z0, s, tau, side); % vel @ test pts
      ue0 = stokesletmatrix(struct('x',z0),y_force)*pt_force; % exact soln
      ue0 = ue0(1:end/2) + 1i*ue0(end/2+1:end); % make complx rep @ test pts
      du0 = ue0-u0; uconst = du0(1); urot = imag(du0(2)-du0(1)); % get 3 params
      uc(ii) = uc(ii) + uconst + 1i*urot*p.x; % correct the evaluated vel field
    end
elseif lp == 'd'
    if side == 'e'
        uc(ii) = StokesDcloseeval(p.x, s, tau, side) + StokesScloseeval(p.x, s, tau, side);
    else
        uc(ii) = StokesDcloseeval(p.x, s, tau, side); % + StokesScloseeval(p.x, s, tau, side);  % for testing alpert errs
    end
end

% plot layer potentials
if v==2
    figure; streamslice(gx,gy,real(uc),imag(uc)); axis equal tight;
    hold on; t = (1:N)'/N*2*pi; plot(real(s.Z(t)),imag(s.Z(t)),'r');
    if lp == 's'
        title('Stokes SLP')
    else
        title('Stokes DLP')
    end
end

% 2.4 compare the SLP soln with the exact soln
if v==1
    figure; imagesc(gx,gy,log10(abs(uc-fhom)));colormap(jet(256));
    colorbar; axis equal tight;
%     title('log_{10} error in u'); 
end
fprintf('#pts: %d',N)
fprintf(', Max error: %.2g \n',max(abs(uc(ii) - fhom(ii))))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function A = stokesletmatrix(t,s) % stokeslet matrix
% t = target seg (x,nx cols), s = src seg
% No jump included on self-interaction.
N = numel(s.x); M = numel(t.x);
d = repmat(t.x, [1 N]) - repmat(s.x.', [M 1]);    % C-# displacements mat
r = abs(d);    % dist matrix R^{MxN}

log_part = kron(eye(2),-log(r));

d1 = real(d)./r;
d2 = imag(d)./r;

cross_part = [d1.^2, d1.*d2; d1.*d2, d2.^2];

A = (log_part + cross_part)/4/pi;

function A = stokesletpmatrix(t,s) % normal derive of stokeslet matrix
% t = target seg (x,nx cols), s = src seg
N = numel(s.x); M = numel(t.x);
d = repmat(t.x, [1 N]) - repmat(s.x.', [M 1]);    % C-# displacements mat
ny = repmat(t.nx, [1 N]);      % identical rows given by src normals
r2 = d.*conj(d);    % dist^2 matrix R^{MxN}

dot_part = -real(ny./d)./r2/pi;
dot_part = repmat(dot_part,2,2);

d1 = real(d);
d2 = imag(d);
cross_part = [d1.^2, d1.*d2; d1.*d2, d2.^2];
A = dot_part.*cross_part;


function A = DLPmatrix(t,s) % double-layer kernel matrix
% t = target seg (x,nx cols), s = src seg
% No jump included on self-interaction.
N = numel(s.x); M = numel(t.x);
d = repmat(t.x, [1 N]) - repmat(s.x.', [M 1]);    % C-# displacements mat
ny = repmat(s.nx.', [M 1]);      % identical rows given by src normals
r2 = d.*conj(d);    % dist^2 matrix R^{MxN}

dot_part = real(ny./d)./r2/pi;
dot_part = repmat(dot_part,2,2);

d1 = real(d);
d2 = imag(d);
cross_part = [d1.^2, d1.*d2; d1.*d2, d2.^2];
A = dot_part.*cross_part;

if numel(s.x)==numel(t.x) && max(abs(s.x-t.x))<1e-14
    t1 = real(s.tang); t2 = imag(s.tang);
    A(diagind(A)) = [t1.^2; t2.^2].*repmat(-s.cur,2,1)/2/pi;
    A = A(:,[1+end/2:end,1:end/2]);
    A(diagind(A)) = [t1.*t2; t2.*t1].*repmat(-s.cur,2,1)/2/pi;
    A = A(:,[1+end/2:end,1:end/2]);
end           % self? diagonal term for Stokes

A = A.*repmat(s.w(:)', [2*M 2]);

function A = SLPmatrixp(t,s) % normal deriv of single-layer kernel self matrix
% t = target seg (x,nx cols), s = src seg
% No jump included on self-interaction.
N = numel(s.x); M = numel(t.x);
d = repmat(t.x, [1 N]) - repmat(s.x.', [M 1]);    % C-# displacements mat
ny = repmat(t.nx, [1 N]);      % identical rows given by src normals
r2 = d.*conj(d);    % dist^2 matrix R^{MxN}

dot_part = -real(ny./d)./r2/pi;
dot_part = repmat(dot_part,2,2);

d1 = real(d);
d2 = imag(d);
cross_part = [d1.^2, d1.*d2; d1.*d2, d2.^2];
A = dot_part.*cross_part;

if numel(s.x)==numel(t.x) && max(abs(s.x-t.x))<1e-14
    t1 = real(s.tang); t2 = imag(s.tang);
    A(diagind(A)) = [t1.^2; t2.^2].*repmat(-s.cur,2,1)/2/pi;
    A = A(:,[1+end/2:end,1:end/2]);
    A(diagind(A)) = [t1.*t2; t2.*t1].*repmat(-s.cur,2,1)/2/pi;
    A = A(:,[1+end/2:end,1:end/2]);
end           % self? diagonal term for Stokes

A = A.*repmat(s.w(:)', [2*M 2]);

function i = diagind(A)
% function i = diagind(A)
%
% return diagonal indices of a square matrix, useful for changing a diagonal
% in O(N) effort, rather than O(N^2) if add a matrix to A using matlab diag()
%
% barnett 2/6/08
N = size(A,1);
if size(A,2)~=N
  disp('input must be square!');
end
i = sub2ind(size(A), 1:N, 1:N);
