function [x_out, w_out] = lglnodes(n, options)
%LGLNODES Compute Legendre-Gauss-Lobatto nodes and weights on a hypercube.
%
%   [x, w] = lglnodes(n) returns 1D nodes and weights.
%   [x, w] = lglnodes(n, "Dimension", D) returns nodes and weights 
%   for a D-dimensional hypercube via tensor products.

arguments
    n (1,1) {mustBeInteger, mustBeGreaterThanOrEqual(n,2)}
    options.Dimension (1,1) {mustBeInteger, mustBePositive} = 1
end

% --- 1. Compute 1D Base Nodes and Weights ---
N = n - 1;
x1D = cos(pi*(0:N)'/N);
xold = 2*ones(size(x1D));
P = zeros(n, n);

% Newton iteration for roots of derivative of Legendre polynomial of degree N
while max(abs(x1D - xold)) > eps
    xold = x1D;
    P(:,1) = 1;
    P(:,2) = x1D;

    for k = 2:N
        P(:,k+1) = ((2*k - 1).*x1D.*P(:,k) - (k - 1).*P(:,k-1))/k;
    end

    x1D = xold - (x1D.*P(:,N+1) - P(:,N))./(n*P(:,N+1));
end

% Return ascending 1D nodes from -1 to 1
x1D = flipud(x1D);

if nargout > 1
    w1D = 2./(N*n*(P(:,N+1).^2));
    w1D = flipud(w1D);
end

% --- 2. Extend to N-Dimensions ---
D = options.Dimension;

if D == 1
    x_out = x1D;
    if nargout > 1
        w_out = w1D;
    end
else
    % Create N-dimensional grid dynamically using cell arrays
    grids = cell(1, D);
    % ndgrid replicates the 1D nodes across all D dimensions
    [grids{:}] = ndgrid(x1D);
    
    % Flatten grids into an (n^D) x D coordinate matrix
    x_out = zeros(n^D, D);
    for d = 1:D
        x_out(:, d) = grids{d}(:);
    end
    
    % Compute tensor product weights if requested
    if nargout > 1
        w_grids = cell(1, D);
        [w_grids{:}] = ndgrid(w1D);
        
        w_out = ones(n^D, 1);
        for d = 1:D
            w_out = w_out .* w_grids{d}(:);
        end
    end
end
end