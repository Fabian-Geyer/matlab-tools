function expr_in_x = xbar_to_x(expr_in_xbar, x_bar_vars, x_vars, x_bar_star, S)
% Transforms a symbolic expression from original (x_bar) to new (x) coordinates.
    if isvector(S), S = diag(S); end
    x_bar_expression = inv(S) * x_vars + x_bar_star(:);
    expr_in_x = subs(expr_in_xbar, x_bar_vars, x_bar_expression);
end
