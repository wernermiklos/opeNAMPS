function out = cell_cross_level1(a,b)
%IRREPSET_CROSS Returns the cartesian product of a, b.
%   a and b are cell arrays whose elements are vectors.
%   out is a cell array with elements out{1} = [a{1},b{1}], out{2} =
%   [a{1},b{2}], ...
    if ~iscell(a) || ~iscell(b)
        error('a and b must be cell arrays!')
    end
    out = cell(1,length(a)*length(b));
    index = 1;
    for i = 1:length(a)
        for j = 1:length(b)
            out{index} = [a{i},b{j}];
            index = index + 1;
        end
    end
end

