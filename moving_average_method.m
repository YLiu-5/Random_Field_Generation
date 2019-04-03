% Subsection 2.6 % 
% 3D Moving Average simulation %
Nx = 10;
Ny = 10;
Nz = 10;

Lx = 1;
Ly = 1;
Lz = 1;

dx = Lx / Nx;
dy = Ly / Ny;
dz = Lz / Nz;

% Correlation length
clx = 1;
cly = 1;
clz = 1;

k_avg = 1;
V_dp = 0.5;

x=linspace(-Lx/2.0,Lx/2.0,Nx);
y=linspace(-Lx/2.0,Lx/2.0,Ny);
z=linspace(-Lx/2.0,Lx/2.0,Nz);
[X,Y,Z]=ndgrid(x,y,z);
s=-log(1-V_dp); % standard deviation
mu=log(k_avg)-s*s/2.0; % mean for a log-random field
Z=s*randn(Nx,Ny,Nz); % normal distribution
F = exp(-(X.^2/(clx*clx/2.0)+Y.^2/(cly*cly/2.0)+Z.^2/(clz*clz/2.0))); 
% Gaussian filter kernel
f = 2.0/sqrt(pi)*Lx/(Nx*Ny*Nz)^(1/3)/(clx*cly*clz)^(1/3).*ifftn(fft2(Z).*fftn(F)); 
K=exp(mu+real(f)); % permeability field
slice(K, [1, Ny], [1, Nx], [1, Nz]);
xlabel('Cell centers [y vlaues]');
ylabel('Cell centers [x vlaues]');
zlabel('Cell centers [z vlaues]');
axis tight
colormap(jet)
colorbar