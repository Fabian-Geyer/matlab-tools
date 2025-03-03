%% One input - one output: Allowed
x = casadi.SX.sym('x');
f = sin(x);
fFun = casadi.Function('f', {x}, {f},{'x'},{'f(x)'});

x0 = 1;

% simulate dynamics
[t, y] = mlt.dynamics.ode45(fFun, x0);


%% One 1x2 input - one output: Allowed
x = casadi.SX.sym('x',2,1);
f = [x(1)-x(2)^2; 
    sin(x(2))];
fFun = casadi.Function('f', {x}, {f},{'x'},{'f(x)'});

x0 = [1;2];

% simulate dynamics
[t, y] = mlt.dynamics.ode45(fFun, x0);


%% Three 1x1 inputs - one 3x1 output: Allowed
x1 = casadi.SX.sym('x1');
x2 = casadi.SX.sym('x2');
x3 = casadi.SX.sym('x3');

f = [x1;x1+x2^2;sqrt(x3)];
fFun = casadi.Function('f', {x1,x2,x3}, {f}, {'x1','x2','x3'}, {'f(x)'});

% simulate dynamics
[t,y] = mlt.dynamics.ode45(fFun, [0,2,3]);

%% One 1x1 input, one 1x2 input - one 3x1 output
x12 = casadi.SX.sym('x',1,2);
x3 = casadi.SX.sym('x3');

f = [x12(1);x12(1)+x12(2)^2;sqrt(x3)];
fFun = casadi.Function('f', {x12,x3}, {f}, {'x','x3'}, {'f(x)'});

% simulate dynamics
[t,y] = mlt.dynamics.ode45(fFun, [0,1,2]);