function [t, X] = euler_maruyama(drift_term, diff_term, x0, T, N, W)
    % EULER_MARUYAMA Simulates an SDE using the Euler-Maruyama method.
    %
    % INPUTS:
    %   drift_term - Drift: CasADi Function f(x), function handle, or constant.
    %   diff_term  - Diffusion: CasADi Function g(x), function handle, or constant.
    %   x0         - Initial condition (numeric or CasADi SX/MX variable).
    %   T          - Total simulation time horizon.
    %   N          - Number of time steps.
    %   W          - (Optional) Standard normal random variables (nw x N).
    %
    % OUTPUTS:
    %   t          - Time vector (1 x N+1)
    %   X          - State trajectory (nx x N+1).

    import casadi.*
    
    dt = T / N;
    t = linspace(0, T, N+1);
    nx = size(x0, 1);
    
    % Evaluate diffusion term at x0 to find the number of noise dimensions (nw)
    g_eval = evaluate_term(diff_term, x0);
    
    if nargin < 6 || isempty(W)
        nw = size(g_eval, 2);
        W = randn(nw, N); 
    else
        nw = size(W, 1);
    end
    
    % Check if we are building a symbolic graph or doing purely numerical simulation
    is_symbolic = isa(x0, 'casadi.MX') || isa(x0, 'casadi.SX') || ...
                  isa(W, 'casadi.MX') || isa(W, 'casadi.SX') || ...
                  isa(g_eval, 'casadi.MX') || isa(g_eval, 'casadi.SX');
    
    if is_symbolic
        % SYMBOLIC UNROLLING
        X_cell = cell(1, N+1);
        X_cell{1} = x0;
        x_curr = x0;
        
        for k = 1:N
            dW_k = W(:, k) * sqrt(dt);
            
            f_val = evaluate_term(drift_term, x_curr);
            g_val = evaluate_term(diff_term, x_curr);
            
            x_next = x_curr + f_val * dt + g_val * dW_k;
            
            X_cell{k+1} = x_next;
            x_curr = x_next;
        end
        
        X = horzcat(X_cell{:});
        
    else
        % NUMERIC EVALUATION
        X = zeros(nx, N+1);
        X(:, 1) = full(x0);
        x_curr = full(x0);
        
        for k = 1:N
            dW_k = W(:, k) * sqrt(dt);
            
            f_val = full(evaluate_term(drift_term, x_curr));
            g_val = full(evaluate_term(diff_term, x_curr));
            
            x_curr = x_curr + f_val * dt + g_val * dW_k;
            X(:, k+1) = x_curr;
        end
    end
end

% --- Helper Function to Handle Constants vs. Functions ---
function val = evaluate_term(term, x)
    % Checks if the term is a callable function or a static constant
    if isa(term, 'casadi.Function') || isa(term, 'function_handle')
        val = term(x);
    else
        val = term; % Return the constant (numeric matrix or symbolic variable)
    end
end