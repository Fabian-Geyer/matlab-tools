function x = cell2vec(C)
%CELL2VEC Flattens all elements in a cell array into a single column vector.
%   Useful for handling CasADi SX/MX symbolic expressions or general cell arrays.
%
%   Usage:
%       x = mlt.utils.cell2vec(f.sx_in);
%
%   Inputs:
%       C - Cell array containing vectors or matrices
%
%   Outputs:
%       x - Column vector containing all concatenated elements

    arguments
        C cell
    end

    % Apply vec to each element and store in a new cell array
    flattenedCells = cellfun(@vec, C, 'UniformOutput', false);

    % Concatenate all elements into a single column vector
    x = vertcat(flattenedCells{:});
end
