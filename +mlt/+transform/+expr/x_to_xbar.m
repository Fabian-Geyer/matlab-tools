function expr_in_xbar = x_to_xbar(expr_in_x, x_vars, x_bar_vars, x_bar_star, S)
% Transforms a symbolic expression from new (x) to original (x_bar) coordinates.
    if isvector(S), S = diag(S); end
    x_expression = S * (x_bar_vars - x_bar_star(:));
    expr_in_xbar = subs(expr_in_x, x_vars, x_expression);
end
