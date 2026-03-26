function [x, w] = lglnodes(n)
%LGLNODES Compute Legendre-Gauss-Lobatto nodes (and optional weights) on [-1, 1].

arguments
    n (1,1) {mustBeInteger, mustBeGreaterThanOrEqual(n,2)}
end

N = n - 1;
x = cos(pi*(0:N)'/N);
xold = 2*ones(size(x));
P = zeros(n, n);

% Newton iteration for roots of derivative of Legendre polynomial of degree N
while max(abs(x - xold)) > eps
    xold = x;
    P(:,1) = 1;
    P(:,2) = x;

    for k = 2:N
        P(:,k+1) = ((2*k - 1).*x.*P(:,k) - (k - 1).*P(:,k-1))/k;
    end

    x = xold - (x.*P(:,N+1) - P(:,N))./(n*P(:,N+1));
end

% Return ascending nodes from -1 to 1
x = flipud(x);

if nargout > 1
    w = 2./(N*n*(P(:,N+1).^2));
    w = flipud(w);
end
end