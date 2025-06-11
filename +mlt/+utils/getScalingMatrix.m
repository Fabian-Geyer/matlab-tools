function S = getScalingMatrix(fun, domX, opts)
%GETSCALINGMATRIX Constructs a diagonal scaling matrix for a CasADi function
%   S = getScalingMatrix(fun, domX, opts) evaluates the CasADi function 'fun'
%   over a hypergrid defined by 'domX' and computes the maximum absolute value
%   for each output component to construct a positive diagonal scaling matrix.

arguments
    fun (1,1) {mustBeA(fun, 'casadi.Function')}
    domX (:,2) double % Each row defines [min, max] for one input
    opts.dimGrid (1,1) {mustBeInteger, mustBePositive} = 5;
end

nInputs = size(domX, 1);
gridVectors = cell(1, nInputs);

% Create linspace vectors for each input variable
for i = 1:nInputs
    gridVectors{i} = linspace(domX(i, 1), domX(i, 2), opts.dimGrid);
end

% Create full grid using ndgrid and reshape to list of points
[gridCells{1:nInputs}] = ndgrid(gridVectors{:});
nGridPoints = numel(gridCells{1});
X = zeros(nInputs, nGridPoints);
for i = 1:nInputs
    X(i, :) = reshape(gridCells{i}, 1, []);
end

% Evaluate function at all points
fun = mlt.utils.toVecInFun(fun);
Y = fun(X); % Evaluate CasADi function
Y = full(Y); % Convert to numeric matrix

% Find max absolute value for each output
maxAbsVals = max(abs(Y), [], 2);

% Construct diagonal scaling matrix
S = diag(maxAbsVals);
end
