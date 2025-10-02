function f_x_bar = x_to_xbar(f_x, x_vars, x_bar_vars, x_bar_star, S)
% Transforms system dynamics from new (x_dot) to original (x_bar_dot) coordinates.
    if isvector(S), S = diag(S); end
    f_in_orig_coords = mlt.transform.expr.x_to_xbar(f_x, x_vars, ...
        x_bar_vars, x_bar_star, S)
    f_x_bar = inv(S) * f_in_orig_coords;
end