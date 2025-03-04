function [vecFun] = toVecInFun(casadiFun)
% Takes a CasADi SX Function with n (1x1) input and returns a CasADi SX
% Function with one vector input of size (nx1)
arguments
    casadiFun (1,1) {mustBeA(casadiFun, 'casadi.Function')}
end

if ~(numel(casadiFun.sx_in)>1)
    warndlg('There is already only one input variable');
    vecFun = casadiFun;
end

% get all inputs
xin = casadiFun.sx_in;

% check if all inputs are 1x1
all1x1 = all(cellfun(@(x) isequal(size(x), [1, 1]), xin));

if ~all1x1
    error('All function inputs must be scalars')
end

% nx1 casadi.SX
xin = mlt.utils.cell2vec(xin);

% Define new function
cellIn = num2cell(xin);
vecFun = casadi.Function(casadiFun.name, {xin}, {casadiFun(cellIn{:})}, ...
    {'x'}, casadiFun.name_out);

end