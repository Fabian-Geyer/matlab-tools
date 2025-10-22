%% Transform point
px = [1;1];
S = diag([1;2]);
xbar_star = [2;2];

% Coordinates should be 0 in x-system:
px0 = mlt.transform.point.xbar_to_x(xbar_star, xbar_star, S)

% (1,1) should be (-1,-2) in x coordinates
px1 = mlt.transform.point.xbar_to_x([1;1], xbar_star, S)

% (0,0) should be xbar_star=(2;2) in bar coordinates
pxbar0 = mlt.transform.point.x_to_xbar([0;0], xbar_star, S)

% (1,1) should be (3,2.5) in bar coordinates
pxbar1 = mlt.transform.point.x_to_xbar([1;1], xbar_star, S)


%% Transform expression
x = casos.PD('x');
xbar = casos.PD('xb');
xbar_star = 1;
S = 2;

% transform xbar to x
expr_x = mlt.transform.expr.xbar_to_x(xbar^2, xbar, x, xbar_star, S)

% transform back
mlt.transform.expr.x_to_xbar(expr_x, x, xbar, xbar_star, S)

%% Transform system
x = casos.PD('x');
xbar = casos.PD('xb');
xbar_star = 1;
S = 2;
xbar_dot = xbar^2

xdot = mlt.transform.system.xbar_to_x(xbar_dot, xbar, x, xbar_star, S)

xbar_dot = mlt.transform.system.x_to_xbar(xdot, x, xbar, xbar_star, S)