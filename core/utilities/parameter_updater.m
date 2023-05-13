function parameters = parameter_updater(parameters,argcell)
%PARAMETER_UPDATER Summary of this function goes here
%   Detailed explanation goes here
for i = 1:2:length(argcell)
    try
        parameters.(argcell{i}) = argcell{i+1};
    catch
        disp(argcell);
        error('varargin parameters should be <name1> val1 <name2> val2, ...');
    end
end

