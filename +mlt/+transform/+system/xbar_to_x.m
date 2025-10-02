function f_x = xbar_to_x(f_x_bar, x_bar_vars, x_vars, x_bar_star, S)
% Transforms system dynamics from original (x_bar_dot) to new (x_dot) coordinates.
    if isvector(S), S = diag(S); end
    f_in_new_coords = mlt.transform.expr.xbar_to_x(f_x_bar, ...
        x_bar_vars, x_vars, x_bar_star, S);
    f_x = S * f_in_new_coords;
end