% Tests for mlt.sde.generator (requires CaSoS on the path).
% Reference: full Hessian, LV = grad(V) f + 1/2 sum(sum((g g') .* Hess(V))).
iz  = @(p) is_zero(remove_coeffs(casos.PD(p), 1e-10));   % zero up to round-off
ref = @(V, f, g, x) nabla(V, x)*f + sum(sum((g*g') .* nabla(nabla(V, x), x)))/2;

%% Scalar GBM: LV = (s^2 - 2a) x^2
a = 2;  s = 0.7;
x = casos.Indeterminates('x');
LV = mlt.sde.generator(x^2, -a*x, s*x, x);
assert(iz(LV - (s^2 - 2*a)*x^2), 'GBM');

%% Additive noise in one channel (Duffing): LV = grad(V) f + s^2/2 d^2V/dx2^2
x = casos.Indeterminates('x', 2);
f = [x(2); -x(1) - x(1)^3 - 0.5*x(2)];
V = x(1)^2 + x(1)*x(2) + x(2)^2 + x(1)^4/2;
[LV, LVd, LVi] = mlt.sde.generator(V, f, [0; s], x);
assert(iz(LVi - s^2), 'Duffing Ito part');
assert(iz(LV - LVd - LVi), 'parts');
assert(iz(LV - ref(V, f, [0; s], x)), 'Duffing');

%% Full numeric g (more columns than states), generic quartic V
f = [-x(1) - x(1)^3 + x(2); -x(2) - x(2)^3 - x(1) + x(1)^2];
g = [1 0.3 0; -0.2 2 0.5];
V = x(1)^4 + 2*x(1)^3*x(2) + 3*x(2)^4 + x(1)*x(2) + x(2)^2;
assert(iz(mlt.sde.generator(V, f, g, x) - ref(V, f, g, x)), 'numeric g');

%% Constant casos.PD g
gp = casos.PD(1).eye(2).*s;
assert(iz(mlt.sde.generator(V, f, gp, x) - ref(V, f, s*eye(2), x)), 'PD constant g');

%% State-dependent g (multiplicative + additive)
g = [x(1) 0.5; x(1)*x(2) 1 + x(2)^2];
assert(iz(mlt.sde.generator(V, f, g, x) - ref(V, f, g, x)), 'polynomial g');

%% Symbolic coefficients (casos.PS), numeric and polynomial g
c = casos.PS.sym('c', [3 1]);
Vs = c(1)*x(1)^2 + c(2)*x(1)*x(2) + c(3)*x(2)^4;
for gg = {[1 0.3; -0.2 2], [x(1) 0.5; 0 1 + x(2)]}
    D = casos.Function('D', {c}, {mlt.sde.generator(Vs, f, gg{1}, x) - ref(Vs, f, gg{1}, x)});
    assert(iz(D([1; -0.4; 2])), 'symbolic coefficients');
end

disp('test_sde_generator: all passed');
