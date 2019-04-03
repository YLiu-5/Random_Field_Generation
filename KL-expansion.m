% Subsection 2.5 %
% 3D KL-expansion simulation %
% The parameter setting and covariance matrix construction are the same as code above (Matrix decomposition method, Subsection 2.4)

% Here we apply the SVD method,
%% SVD,
% Parameter Setting %
nx = 10;
ny = 10;
nz = 10;

Lx = 1;
Ly = 1;
Lz = 1;

dx = Lx / nx;
dy = Ly / ny;
dz = Lz / nz;
covMat=zeros(nx*ny);

% Matern covariance, with parameters 
nu=1;
lambda=1;
sigma2=1;

% Measure distance between ii and jj, construct covariance matrix
for ii=1:nx*ny*nz
    for jj=1:nx*ny*nz
        pageii = floor(ii/(nx*ny));
        premi   = rem(ii,(nx*ny));
        rowii  = floor(premi/nx);
        colii  = rem(premi,nx);
        
        pagejj = floor(jj/(nx*ny));
        premj   = rem(jj,(nx*ny));
        rowjj  = floor(premj/nx);
        coljj  = rem(premj,nx);
        
        z=abs((pagejj-pageii))*dz;
        y=abs((rowjj-rowii))*dy;
        x=abs((coljj-colii))*dx;
        d=(z^2+y^2+x^2)^(1/2);
        if (d==0)
            covMat(ii,jj)=sigma2;
        else
        covMat(ii,jj)=sigma2*2^(1-nu)/gamma(nu)*(sqrt(2*nu)*d/lambda)^(nu)*besselk(nu,sqrt(2*nu)*d/lambda);
        end
    end
end

Nkl = 300;
[eigenVec,eigenVal]=eigs(covMat,Nkl);
eigenVal=diag(eigenVal);
xi = randn(Nkl,1);
K = eigenVec * (sqrt(eigenVal).*xi);

Kr = reshape(K,[nx ny nz]);
Kre = exp(Kr);

slice(Kre, [1, ny], [1, nx], [1, nz]);
xlabel('Cell centers [y vlaues]');
ylabel('Cell centers [x vlaues]');
zlabel('Cell centers [z vlaues]');
axis tight
colormap(jet)
colorbar