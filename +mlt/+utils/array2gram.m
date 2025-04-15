function [Q] = array2gram(array,dim)
%ARRAY2GRAM returns cell array of gram matrices from array and dim vector
%containing dimensions of each gram matrix (number of rows)
arguments
    array (:,1) {mustBeNumeric}
    dim (:,1) {mustBeInteger, mustBePositive}
end

% Split array into gram matrix arrays (save as cell array)
Q = mat2cell(array, dim.^2);

% Turn each array into nxn matrix
Q = cellfun(@(x) reshape(x, sqrt(numel(x)), []), Q, UniformOutput=false);
end

