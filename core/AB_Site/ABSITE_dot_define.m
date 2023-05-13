function obj = ABSITE_dot_define(obj_in, name1, name2, newname, varargin)
            % Defines a new operator from operators name1 and name2.
            obj = obj_in;
            parameters = struct();
            parameters.SYMMETRIES = {};
            parameters = parameter_updater(parameters,varargin);

            obj.fermionicities(newname) = obj.fermionicities(name1)*obj.fermionicities(name2);
            obj.operators(newname) = NTdot(obj.operators(name1),obj.operators(name2),{'tau'},{'tau~'});
            if isKey(obj.operator_irreps, name1) && isKey(obj.operator_irreps, name2)    %op irreps are set ==> we can define the op irrep of the product
                                                                                         %if SYMMETRIES are provided (varargin)
                if isempty(parameters.SYMMETRIES)
                    warning(['Opirrep is defined both for "', name1, '" and "', name2, '". Product operators opirrep can only', ...
                             ' be determined if "SYMMETRIES" varargin argument is also provided.']);
                else
                    prodrep = SUPER_fusion_rule(parameters.SYMMETRIES,obj.operator_irreps(name1),obj.operator_irreps(name2));
                    obj.operator_irreps(newname) = prodrep{1}{1};
                end
            end
end

