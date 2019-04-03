% 2D Turning Bands simulation %

% Parameter Setting %
nx = 10;
ny = 10;

dx = 0.1;
dy = 0.1;

sigma2 = 1;  % Variance
lambda = 1;  % Correlation length

% Center of each block
coord0 = [ones(nx,1)*0 ((0:ny-1)+0.5)'*dy];
coord = kron(coord0,ones(nx,1)) + kron(ones(nx,1),((0:nx-1)+0.5)'*[dx 0]);
coord = coord - [nx*dx/2 ny*dy/2];

% 2D lines, orientations determined by van der Corput, or Halton sequence

NL = 16;

% 2 dimensional halton set
H = haltonset(2);
h = net(H,NL);

% lines orientation
lines = [cos(2*pi*h(:,1)) sin(2*pi*h(:,1))];

% Projection

P = coord*(lines)';  %nx*ny by NL matrix, each row is the coordinates of a point projected to NL lines

% Transform 2D covariance into radical spectral density function and then
% the corresponding 1D spectral density
% Assume exponential covariance C(h) = exp ( -h/lambda)


% Generate line process by FFT

% Generate dW

deltak = 0.5;
M = 100;

phi = rand(NL,M);

j = (0:M-1);
omega = j * deltak;
S = (sigma2/2.*omega*lambda).^2/(1+(omega.*lambda).^2).^(3/2);

dW = zeros(NL,M);
for i = 1:NL
    dW(i,:) = (4*S*deltak)^(1/2).*cos(2*pi*phi(i,:)) + 1i*(4*S*deltak)^(1/2).*sin(2*pi*phi(i,:));
end

% Use fft to generate X, which is complex line process
X = fft(dW');
% Extract the real part
Z = real(X);

deltazeta = 1/(M*deltak);
index = zeros(nx*ny,NL);
for i = 1:nx*ny
    index(i,:) = floor(P(i,:)/deltazeta);
end

index = index + M/2; 

% for each point
K = zeros(nx*ny,1);
for i = 1:nx*ny
    idx = index(i,:);
    sum = 0;
    for j = 1:NL
        sum = sum + Z(idx(j),j);
    end
    K(i) = sum/sqrt(NL);
end

Kr = reshape(K,[nx,ny]);

contourf(Kr)