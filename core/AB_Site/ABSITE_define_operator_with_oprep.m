function obj = ABSITE_define_operator_with_oprep(obj, ...
                                                 name, ...                  % char array, the name of the operator (like 'fdag')
                                                 operator, ...              % operator matrix (NAtensor with 'tau' and 'tau~' legs)
                                                 opirrep, ...               % Representation of the operator
                                                 varargin)              
             
    parameters = struct();
    parameters.FERMIONICITY = +1;
    parameters.SYMMETRIES = {};
    parameters = parameter_updater(parameters,varargin);
    
    obj.operators(name) = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},obj.no_of_symmetries);
    obj.operators(name) = NTadd(obj.operators(name), operator);
    obj.operator_irreps(name) = opirrep;
    obj.fermionicities(name) = parameters.FERMIONICITY;

    %if SYMMETRIES is passed, we can make a test
    if ~isempty(parameters.SYMMETRIES)
        irreps = NTget_active_irrep_values(obj.operators(name),{{'tau',1},{'tau~',1}});
        for bID = 1:length(irreps)
            if ~SUPER_selection_rule(parameters.SYMMETRIES,irreps{bID}{1},opirrep,irreps{bID}{2})
                warning(['Selection rule is not satisfied. Unexpected block: ', num2str([irreps{bID}{1}, irreps{bID}{2}]), ', opirrep:', num2str(opirrep)]);
            end
        end
    end
end

