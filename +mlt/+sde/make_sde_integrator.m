function F = make_sde_integrator(nx, f, g, dt, opts)
%MAKE_SDE_INTEGRATOR  Build a reusable one-step SDE integrator (CasADi).
%
%   F = make_sde_integrator(nx, f, g, dt)
%   F = make_sde_integrator(nx, f, g, dt, Name=Value, ...)
%
%   Integrates   dX = f(X) dt + g(X) dW,    dW ~ N(0, dt * I_nw)
%   and returns a CasADi Function:
%
%       x_next = F(x, dW)              x : nx x 1 ,  dW : nw x 1
%
%   Compose it:
%       F.mapaccum(N)   -> roll a full path:  (x0, DW[nw x N]) -> X[nx x N]
%       F.map(P)        -> vectorize over P independent evaluations
%
%   Required (positional)
%     nx    state dimension           (positive integer)
%     f     drift,     @(x) -> nx x 1
%     g     diffusion, @(x) -> nx x nw
%     dt    time step                 (positive scalar)
%
%   Optional (Name=Value)
%     nw       Wiener dimension. Default: inferred from columns of g(x).
%     scheme   "euler-maruyama" (default) | "milstein" | "heun"
%     project  post-step map @(x) -> nx x 1. Default: identity.
%              e.g. @(q) q/sqrt(q'*q) for unit quaternions.
%     name     CasADi Function name.  Default: "F_sde".
%
%   Scheme notes
%     euler-maruyama : strong order 0.5 (1.0 for additive noise). Robust default.
%     milstein       : strong order 1.0; uses CasADi jacobian(g,x). Scalar noise
%                      (nw = 1) only. Ito interpretation.
%     heun           : predictor-corrector, derivative-free, strong order 1.0 for
%                      scalar/commutative noise. STRATONOVICH interpretation --
%                      if your model is in Ito form, feed the corrected drift
%                      f_Strat = f_Ito - 0.5 * sum_j (d g_j / dx) * g_j.

    arguments
        nx  (1,1) double        {mustBeInteger, mustBePositive}
        f   (1,1) function_handle
        g   (1,1) function_handle
        dt  (1,1) double        {mustBePositive}
        opts.nw       double    {mustBeScalarOrEmpty, mustBeInteger, mustBePositive} = []
        opts.scheme   (1,1) string {mustBeMember(opts.scheme, ...
                          ["euler-maruyama","milstein","heun"])} = "euler-maruyama"
        opts.project  (1,1) function_handle = @(x) x
        opts.name     (1,1) string = "F_sde"
    end

    % --- symbolic state and user dynamics -------------------------------
    x  = casadi.SX.sym('x', nx);
    f0 = casadi.SX(f(x));          % cast so constant (additive) g/f are fine too
    g0 = casadi.SX(g(x));

    % --- structural checks on what f and g actually returned ------------
    assert(size(f0,1) == nx && size(f0,2) == 1, ...
        'f(x) must return an %dx1 column vector; got %dx%d.', ...
        nx, size(f0,1), size(f0,2));
    assert(size(g0,1) == nx, ...
        'g(x) must return %d rows; got %d.', nx, size(g0,1));

    nw = size(g0, 2);
    if ~isempty(opts.nw)
        assert(opts.nw == nw, ...
            'Declared nw = %d does not match g(x) column count %d.', opts.nw, nw);
    end

    dW = casadi.SX.sym('dW', nw);
    h  = dt;

    % --- assemble one step of the chosen scheme ------------------------
    switch opts.scheme
        case "euler-maruyama"
            x_upd = x + f0*h + g0*dW;

        case "milstein"
            assert(nw == 1, ['milstein is implemented for scalar noise ', ...
                '(nw = 1). For multi-dim noise use "heun" or full Milstein ', ...
                'with Levy areas.']);
            Jg    = casadi.jacobian(g0, x);              % nx x nx
            x_upd = x + f0*h + g0*dW + 0.5*(Jg*g0).*(dW.^2 - h);

        case "heun"
            xb    = x + f0*h + g0*dW;                    % Euler-Maruyama predictor
            fb    = casadi.SX(f(xb));
            gb    = casadi.SX(g(xb));
            x_upd = x + 0.5*(f0 + fb)*h + 0.5*(g0 + gb)*dW;
    end

    x_next = opts.project(x_upd);
    F = casadi.Function(char(opts.name), {x, dW}, {x_next}, {'x','dW'}, {'x_next'});
end