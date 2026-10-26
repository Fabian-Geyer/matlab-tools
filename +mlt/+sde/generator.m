function [LV, LV_drift, LV_ito] = generator(V, f, g, x)
    % GENERATOR Infinitesimal generator of an Ito SDE applied to a polynomial.
    %
    %   For dX = f(X) dt + g(X) dB, B an m-dimensional Brownian motion,
    %
    %       LV = grad(V)' f + 1/2 trace(g' Hess(V) g).
    %
    % INPUTS:
    %   V - Scalar polynomial (casos.PD or casos.PS, may have symbolic coefficients).
    %   f - Drift, n x 1 polynomial.
    %   g - Diffusion, n x m: numeric matrix or polynomial (casos.PD / casos.PS).
    %   x - Indeterminates, n entries (casos.Indeterminates).
    %
    % OUTPUTS:
    %   LV       - Generator applied to V.
    %   LV_drift - Drift part grad(V)' f.
    %   LV_ito   - Ito (second-order) part 1/2 trace(g' Hess(V) g).
    %
    % If g does not depend on x (numeric or zero degree), the Ito part is
    % computed column-wise as 1/2 sum_k g_k' grad(grad(V)' g_k), which needs
    % no Hessian and skips zero noise columns. Otherwise only the Hessian
    % entries with a nonzero entry in g g' are computed (upper triangle).
    %
    % EXAMPLE (GBM, LV = (s^2 - 2a) x^2):
    %   x = casos.Indeterminates('x');
    %   LV = mlt.sde.generator(x^2, -a*x, s*x, x);

    n = numel(x);
    assert(isscalar(V), 'mlt:sde:generator:V', 'V must be scalar.');
    assert(numel(f) == n, 'mlt:sde:generator:f', 'f must have %d entries (one per indeterminate).', n);
    assert(size(g, 1) == n, 'mlt:sde:generator:g', 'g must have %d rows (one per indeterminate).', n);
    f = reshape(f, n, 1);

    dV = nabla(V, x);                                  % 1 x n gradient
    LV_drift = dV*f;

    if isnumeric(g) || is_zerodegree(g)
        % State-independent diffusion: directional second derivatives
        if isnumeric(g), g = full(g); end
        LV_ito = casos.PD(0);
        for k = 1:size(g, 2)
            gk = g(:, k);
            if isnumeric(gk) && ~any(gk), continue; end
            LV_ito = LV_ito + nabla(dV*gk, x)*gk;
        end
        LV_ito = LV_ito/2;
    else
        % State-dependent diffusion: Q = g g', only the needed Hessian entries
        Q = g*g';
        LV_ito = casos.PD(0);
        for i = 1:n
            for j = i:n
                if is_zero(Q(i, j)), continue; end
                w = 1 + (j > i);                        % off-diagonal entries count twice
                LV_ito = LV_ito + w*Q(i, j)*nabla(dV(i), x(j));
            end
        end
        LV_ito = LV_ito/2;
    end

    LV = LV_drift + LV_ito;
end
