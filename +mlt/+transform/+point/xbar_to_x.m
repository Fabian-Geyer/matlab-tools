function p_x = xbar_to_x(p_xbar, x_bar_star, S)
% Maps a numeric point from original (x_bar) to new (x) coordinates.
    if isvector(S), S = diag(S); end
    p_x = S * (p_xbar(:) - x_bar_star(:));
end