function p_xbar = x_to_xbar(p_x, x_bar_star, S)
% Maps a numeric point from new (x) to original (x_bar) coordinates.
    if isvector(S), S = diag(S); end
    p_xbar = inv(S) * p_x(:) + x_bar_star(:);
end