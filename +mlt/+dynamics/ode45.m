function [t,y] = ode45(f, x0, opts)
% Simulate dynamics xdot = f(x) for CasADi function with n inputs and
% outputs. 
% Allowed function inputs: 
%   - cell array with n elements: each [1x1] casadi.SX
%   - cell array with 1 element: [nx1] casadi.SX
% Allowed function output:
%   - cell array with 1 element: [nx1] casadi.SX

arguments
    f (1,1) {mustBeA(f, 'casadi.Function')}
    x0 (:,1) {mustBeNumeric}
    opts.duration__s (1,1) {mustBeNumeric} = 500;
    opts.n_steps (1,1) {mustBeInteger} = 500;
end


% Check function output
if ~(size(f)==[1,1])
    error('dynamics function output must be single vector')
end

% Number of outputs
nout = numel(f.sx_out{1});

% Check inputs
xin = f.sx_in;
if (isscalar(xin)) % input is vector
    % make sure it is a column vector input
    if ~(size(xin{1},2)==1)
        error('CasADi dynamics function input must be column vector')
    end
    xin = vec(xin{1});
else
    all1x1 = all(cellfun(@(x) isequal(size(x), [1, 1]), xin));
    if ~all1x1
        error('Dynamics function input variables have varying sizes')
    end
    xin = mlt.utils.cell2vec(xin); % nx1 casadi.SX

    % Make sure function input is vector
    cellIn = num2cell(xin);
    f = casadi.Function('f', {xin}, {f(cellIn{:})}, {'x'}, {'xdot'});
end



if ~(nout==size(xin))
    error('Number of inputs varies from number of outputs');
end

%% Runge Kutta 45

x = xin;

tf = opts.duration__s;
h = tf/opts.n_steps;

% RK 45
rk = 4; % subintervals
dt = h/rk;
x0sym = casadi.SX.sym('x0',size(x));

k1 = casadi.SX(length(x),rk); 
k2 = casadi.SX(length(x),rk); 
k3 = casadi.SX(length(x),rk); 
k4 = casadi.SX(length(x),rk);
xk = [x0sym casadi.SX(length(x),rk)];

for j=1:rk
    k1(:,j) = f(xk(:,j));
    k2(:,j) = f(xk(:,j) + dt*k1(:,j)/2);
    k3(:,j) = f(xk(:,j) + dt*k2(:,j)/2);
    k4(:,j) = f(xk(:,j) + dt*k3(:,j));
    xk(:,j+1) = xk(:,j) + dt*(k1(:,j) + 2*k2(:,j) + 2*k3(:,j) + k4(:,j))/6;
end

% integrator
xkp1 = casadi.Function('fk', {x0sym}, {xk(:,end)}, {'x0'}, {'xk'});

y = zeros(length(x0),opts.n_steps+1);
y(:,1) = x0;

% Simulate
for k = 1:opts.n_steps
    y(:,k+1) = full(xkp1(y(:,k)));
end

t = linspace(0,opts.duration__s,size(y,2));

end