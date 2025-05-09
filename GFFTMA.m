function [datasim, G, U] = GFFTMA(model, c, nu, param, seed, nbsimul, nx, dx, ny, dy, nz, dz)
% GFFTMA: Function to simulate non-LMC (but symmetrical) models using FFT-MA
%
% Syntax: 
%   datasim = GFFTMA(model, seed, nbsimul, vsiz, nx, dx, ny, dy, nz, dz)
%
% Description:
%   Simulates spatial fields using FFT-MA. The method supports 1D, 2D, or 3D simulations.
%   - 1D: Specify nx and dx.
%   - 2D: Specify nx, dx, ny, and dy.
%   - 3D: Specify all parameters (nx, dx, ny, dy, nz, dz).
%
% Inputs:
%   model   : nvar x nvar cell array of structures defining the covariance models.
%             - Each cell {i, j} contains structures for cross-covariance 
%               between variables i and j.
%             - Isotropic case: [model type, range, shape, sill].
%             - Anisotropic case: [model type, range1, range2, range3, rotx, roty, rotz, shape, sill].
%   c       : nvar x nvar cell array of sill defining the covariance models.
%   nu      : nvar x nvar cell array of shape parameters defining the covariance models (i.e., for Matern covariances).
%   param   : Contains important parameter such as conditioning data.
%             param.HD is a cell nvar by 1 containing location and value
%             for contioning
%   seed    : Random seed for reproducibility.
%   nbsimul : Number of realizations to generate.
%   vsiz    : Vector specifying [Nx, Ny, Nz] for the simulation grid size.
%             If empty, the grid size is computed automatically to avoid aliasing.
%   nx, dx  : Number of grid points and step size along x-axis.
%   ny, dy  : (Optional) Number of grid points and step size along y-axis.
%   nz, dz  : (Optional) Number of grid points and step size along z-axis.
%
% Outputs:
%   datasim : Cell array of size [nbsimul, nvar]. Each cell contains a vector of
%             simulated values for each variable in the grid.
%
% Author:
%   D. Marcotte, May 2015

% Determine simulation dimensionality
if nargin == 8
    cas = 1; % 1D simulation
    ny = 1; nz = 1; dy = 1; dz = 1;
elseif nargin == 10
    cas = 2; % 2D simulation
    nz = 1; dz = 1;
elseif nargin == 12
    cas = 3; % 3D simulation
else
    error('Check input parameters: Incorrect number of inputs.');
end

% Initialize random seed
rng('default');
rng(round(seed));

% Number of variables
nvar = size(model, 1);

% Enlarge grid size to avoir aliasing
Nx=nx*2+1; Ny=ny*2+1; Nz=nz*2+1;

% Compute grid midpoints
Nx2 = floor(Nx / 2);
Ny2 = floor(Ny / 2);
Nz2 = floor(Nz / 2);

% Generate grid coordinates
if cas == 1
    x0 = (-(Nx2 + 1) * dx : dx : Nx2 * dx)';
    x0r = 0;
elseif cas == 2
    x0 = grille2(-(Nx2 + 1) * dx, Nx2 * dx, dx, -(Ny2 + 1) * dy, Ny2 * dy, dy);
    x0r = [0 0];
elseif cas == 3
    x0 = grille3(-(Nx2 + 1) * dx, Nx2 * dx, dx, -(Ny2 + 1) * dy, Ny2 * dy, dy, -(Nz2 + 1) * dz, Nz2 * dz, dz);
    x0r = [0 0 0];
end

% Initialize covariance matrix
G = ComputeS(x0, x0r, model, c, nu, Nx, Ny, Nz);

% Initialize output data
[n, p] = size(G{1, 1});
datasim = zeros(nx*ny*nz, nbsimul, nvar);
U = cell(nbsimul, 1);
% Perform simulations
for is = 1:nbsimul
    GU = cell(nvar, 1);
    for i = 1:nvar
        GU{i} = 0;
    end
    for j = 1:nvar
        U{is}{j} = randn(n, p); % Fourier transform of random numbers
        for i = 1:nvar
            GU{i} = GU{i} + G{i, j} .* fftn(U{is}{j});
        end
    end
    for j = 1:nvar
        datasimTemp = real(ifftn(GU{j})); % Simulated values for variable j
        datasim(:,is,j) = reshape(datasimTemp(1:nx, 1:ny, 1:nz),[],1);
    end
end

% post-conditioning by simple cokriging
if ~isempty(param.HD)
    datasim = postcond(datasim, param, model, c, nu, nx, ny, nz);
end

end

function G = ComputeS(x0, x0r, model, c, nu, Nx, Ny, Nz)
% ComputeS: Function to compute the spectral densities of multivariate covariance matrices
% and perform Cholesky decomposition for spectral simulation.
%
% Syntax:
%   G = ComputeS(x0, x0r, model, Nx, Ny, Nz)
%
% Description:
%   This function calculates the spectral decomposition of multivariate covariance matrices
%   for a given set of input parameters. It uses the covariance model to compute
%   Fourier-transformed covariance matrices and applies Cholesky decomposition
%   at all frequencies for spectral simulation.
%
% Inputs:
%   x0     : [n x d] matrix of input coordinates.
%   x0r    : [1 x d] matrix of the coordinates of the center grid.
%   model  : nvar x nvar cell array of structures defining the covariance models.
%             - Each cell {i, j} contains structures for cross-covariance
%               between variables i and j.
%             - Isotropic case: [model type, range, shape, sill].
%             - Anisotropic case: [model type, range1, range2, range3, rotx, roty, rotz, shape, sill].
%   nn     : Integer, number of realizations.
%   Nx     : Integer, number of grid points along the x-axis (for 1D cases).
%   Ny     : Integer, number of grid points along the y-axis (for 2D cases).
%   Nz     : Integer, number of grid points along the z-axis (for 3D cases).
%
% Outputs:
%   G      : [nvar x nvar] cell array containing the Cholesky-decomposed
%            spectral matrices at all frequencies.
%
% Notes:
%   - Supports 1D, 2D, or 3D cases based on grid dimensions.
%   - Uses FFT to compute the Fourier transform of the covariance matrices.
%   - Handles isotropic and anisotropic covariance models.
%
% Author:
%   D. Marcotte, May 2015.
%   D. Lauzon , November 2024. (Modified from Marcotte, 2015)
%

% Compute the case
cas = size(x0,2);

% Number of variables
nvar = size(model, 1);

% Initialize covariance matrix
S = cell(nvar, nvar);

for i = 1:nvar
    for j = 1:i
        cc = covar(x0, x0r, model{i, j}, c{i, j}, nu{i, j});

        % Handle different cases
        switch cas
            case 1
                cc = fftshift(cc);
            case 2
                cc = fftshift(reshape(cc, Nx + 1, Ny + 1));
            case 3
                cc = fftshift(reshape(cc, Nx + 1, Ny + 1, Nz + 1));
        end

        % Fourier transform of covariances
        S{i, j} = real(fftn(cc));
        S{j, i} = S{i, j}; % Symmetry
    end
end

% Perform Cholesky decomposition
G = chol_dec(S);

end

function G = chol_dec(S)
% chol_dec: Helper function for Cholesky decomposition of spectral matrices.
%
% Syntax:
%   G = chol_dec(S)
%
% Description:
%   Performs Cholesky decomposition of spectral matrices for use in
%   spectral simulation. This function handles matrices stored in a
%   cell array format and processes only the significant frequencies.
%
% Inputs:
%   S      : [nvar x nvar] cell array of spectral covariance matrices.
%            - Each cell contains a matrix of spectral coefficients
%              corresponding to a given variable pair.
%
% Outputs:
%   G      : [nvar x nvar] cell array of decomposed spectral matrices.
%            - Each cell contains a matrix of Cholesky decomposition
%              results at all frequencies.
%
% Notes:
%   - Applies eigenvalue decomposition to ensure non-negative eigenvalues.
%   - Processes only significant frequencies to improve computational efficiency.
%
% Author:
%   D. Marcotte, May 2015.
%   D. Lauzon , November 2024. Vectorized the code.
%
% Perform Cholesky decomposition simultaneously at all frequencies
% S is a cell array (nvar x nvar), each cell containing a matrix of spectrum coefficients

nvar = size(S, 1);
[nx, ny] = size(S{1, 1});
total_points = nx * ny;

% Identify significant indices
id = zeros(nx, ny);
for i = 1:nvar
    id = id | (S{i, i} > max(S{i, i}(:)) * 1e-6); % Logical OR for significant values
end
id = find(id); % Linear indices of significant points
num_points = length(id);

% Preallocate 3D array for significant indices: nvar x nvar x num_points
m = zeros(nvar, nvar, num_points);

% Populate m matrices for significant indices
for i1 = 1:nvar
    for i2 = i1:nvar
        values = reshape(S{i1, i2}, total_points, 1);
        m(i1, i2, :) = values(id);
        if i1 ~= i2
            m(i2, i1, :) = values(id);
        end
    end
end

% Perform eigenvalue decomposition in batch
[V, D] = pageeig(m); % Batch eigenvalue decomposition (nvar x nvar x num_points)

% Construct G using eigen decomposition results
D = max(D, 0); % Ensure non-negative eigenvalues
G_matrices = pagemtimes(V, pagemtimes(sqrt(D), pagetranspose(V))); % G = V * sqrt(D) * V'

% Reassign results back to G
G = cell(nvar, nvar);
for i1 = 1:nvar
    for i2 = 1:nvar
        temp = zeros(total_points, 1);
        temp(id) = squeeze(G_matrices(i1, i2, :));
        G{i1, i2} = reshape(temp, nx, ny);
    end
end
end

function datasim = postcond(datasim, param, model, c, nu, nx, ny, nz)
% postcond: Performs post-conditioning by simple cokriging.
%
% Inputs:
%   datasim - Simulated data matrix (N x nVar).
%   param   - A structure containing conditioning parameters, including param.HD.
%   model   - Covariance model definitions.
%   c       - Covariance model sills.
%   nu      - Covariance model shape parameters.
%   nx, ny, nz - Grid dimensions.
%
% Output:
%   paramStruct - A structure containing all necessary information for post-conditioning.

% Determine simulation dimensionality and grid points
if nz == 1 && ny == 1 && nx > 1
    x0 = (1:nx)';  % 1D grid
    d = 1;
elseif nz == 1 && ny > 1 && nx > 1
    x0 = grille2(1, nx, 1, 1, ny, 1);  % 2D grid
    d = 2;
elseif nz > 1 && ny > 1 && nx > 1
    x0 = grille3(1, nx, 1, 1, ny, 1, 1, nz, 1);  % 3D grid
    d = 3;
else
    error('Check input parameters: Incorrect grid dimensions.');
end

nvar = size(model, 1);  % Number of variables
x0grid = cell(nvar, 1); x = cell(nvar, 1); xgrid = cell(nvar, 1);
idx = cell(nvar, 1); xsim = cell(nvar, 1); x0sim = cell(nvar, 1);

% % Permute datasim to be efficient with a previous code
% datasim = permute(datasim,[1 3 2]);

% Loop over the number of variables
for i = 1:nvar
    x{i, 1} = param.HD{i}(~isnan(param.HD{i}(:, d + 1)), d + 1);
    xgrid{i, 1} = param.HD{i}(~isnan(param.HD{i}(:, d + 1)), 1:2);
    idx{i} = (xgrid{i, 1}(:, 2) - 1) * nx + xgrid{i, 1}(:, 1);
    xsim{i, 1} = datasim(idx{i},:, i);
    x0sim{i, 1} = datasim(:,:, i);
end

% Precompute covariance matrices of observed points
K = ComputeS_Cokri(xgrid, xgrid, model, c, nu, 1);

% Extract precomputed values
nbsimul = size(x0sim{1},2);

% Vectorize inputs
x = cell2mat(x);
xsim = cell2mat(xsim);

% Compute residuals
xres = x - xsim;

% Solve the kriging system (dual form)
bres = K \ xres;

% Iteratively cokriged the residuals by batch
batchSize = 10000; % Define the size of each batch
M = size(x0,1); % Total number of simulated point
numBatches = ceil(M/ batchSize); % Calculate the total number of batches

%Reinitialized datasim
datasim = zeros(M,nbsimul,nvar);
for i = 1 : nvar
    for batchIdx = 1:numBatches
        % Calculate the start and end column indices for the current batch
        startCol = (batchIdx - 1) * batchSize + 1;
        endCol = min(batchIdx * batchSize, M);        
        for j=1:nvar
            if i == j
                x0grid{j}=x0(startCol:endCol,:);
            else
                x0grid{j}=[];
            end
        end
        % Covariance matrices of observed points with grid point
        k = ComputeS_Cokri(xgrid, x0grid, model, c, nu, 0);

        % Compute cokriged residuals
        xres_krig = (bres' * k)';

        % Combine unconditional simulations and cokriged residuals
        datasim(startCol:endCol,:,i) = xres_krig + x0sim{i}(startCol:endCol,:);
    end
end

end

function C = ComputeS_Cokri(x0, x0r, model, c, nu, idx)
% ComputeS: Function to compute the spectral densities of multivariate covariance matrices
% and perform Cholesky decomposition for spectral simulation.
%
% Syntax:
%   G = ComputeS_Cokri(x0, x0r, model)
%
% Description:
%   This function calculates the spectral decomposition of multivariate covariance matrices
%   for a given set of input parameters. It uses the covariance model to compute
%   Fourier-transformed covariance matrices and applies Cholesky decomposition
%   at all frequencies for spectral simulation.
%
% Inputs:
%   x0     : [n x d] matrix of input coordinates.
%   x0r    : [1 x d] matrix of the coordinates of the center grid.
%   model  : nvar x nvar cell array of structures defining the covariance models.
%             - Each cell {i, j} contains structures for cross-covariance
%               between variables i and j.
%             - Isotropic case: [model type, range, shape, sill].
%             - Anisotropic case: [model type, range1, range2, range3, rotx, roty, rotz, shape, sill].
%   nn     : Integer, number of realizations.
%
% Outputs:
%   G      : [nvar x nvar] cell array containing the Cholesky-decomposed
%            spectral matrices at all frequencies.
%
% Notes:
%   - Supports 1D, 2D, or 3D cases based on grid dimensions.
%   - Uses FFT to compute the Fourier transform of the covariance matrices.
%   - Handles isotropic and anisotropic covariance models.
%
% Author:
%   D. Marcotte, May 2015.
%   D. Lauzon , November 2024. (Modified from Marcotte, 2015)
%


% Number of variables
nvar = size(model, 1);

% Initialize covariance matrix
C = cell(nvar, nvar);

% Perform Cholesky decomposition
if idx == 1
    [v , d] = eig(cell2mat(C));
    C = v * diag(max(diag(d), 0)) * v';
    for i = 1:nvar
        for j = 1:i
            C{i, j} = covar(x0{i}, x0r{j}, model{i, j}, c{i, j}, nu{i, j});
            C{j, i} = C{i, j}'; % Symmetry
        end
    end
    C = cell2mat(C);
else
    for i = 1:nvar
        for j = 1:nvar
            if isempty(x0r{j})
                C{i, j}= [];
            else
                C{i, j} = covar(x0{i}, x0r{j}, model{i, j}, c{i, j}, nu{i, j});
            end
        end
    end
    C = cell2mat(C);
end
end




