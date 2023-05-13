function obj = SITE_create(symmetries,sector_list,varargin)
            % Constructor of "Site" objects
            % symmetries:              cell array of Symmetry objects
            % sector_list:             cell array of sectors
            %                   format: {{[irrep labels1], taudim1},
            %                   {[irrep labels2], taudim2}, ...}
            obj.type = 'Site';
            
            parameters = struct();
            % Default parameters
            parameters.FERMIONS = false;
            % Update according to varargin
            parameters = parameter_updater(parameters,varargin);

            obj.fermions = parameters.FERMIONS;

            
            obj.no_of_symmetries = length(symmetries);
            obj.sector_multiplets = containers.Map();
            obj.multiplet_dims = containers.Map();
            for secID = 1:length(sector_list)
                if length(sector_list{secID}{1}) ~= obj.no_of_symmetries
                    error('Irrep labels must contain no_of_symmetries quantum numbers');
                end
                seckey = char(sector_list{secID}{1});
                irrepdim = 1;
                for symID = 1:obj.no_of_symmetries
                    irrepdim = irrepdim * symmetries{symID}.irrep_dimensions(sector_list{secID}{1}(symID));  
                end
                obj.sector_multiplets(seckey) = sector_list{secID}{2};
                obj.multiplet_dims(seckey) = irrepdim;
            end
            obj.operators = containers.Map();
            obj.fermionicities = containers.Map();
            obj.op_is_scalar = containers.Map();
            obj.eta_legs = containers.Map();
            obj.operators('id') = NAtensor({'tau','tau~'},{'o','i'},{[1],[1]},obj.no_of_symmetries);
            for seckey_cell = obj.sector_multiplets.keys
                seckey = seckey_cell{1};
                obj.operators('id') = NTset_block(obj.operators('id'),{{'tau',1}},{double(seckey)},{'tau','tau~'},eye(obj.sector_multiplets(seckey)));
            end
            obj.fermionicities('id') = 1;
            obj.op_is_scalar('id') = true;
            obj.eta_legs('id') = {};
            if parameters.FERMIONS   % We create a dummy phase operator (should be overwritten in general, but works for the auxiliary site of DMRG)
                obj.operators('ph') = NAtensor({'tau','tau~'},{'o','i'},{[1],[1]},obj.no_of_symmetries);
                for seckey_cell = obj.sector_multiplets.keys
                    seckey = seckey_cell{1};
                    obj.operators('ph') = NTset_block(obj.operators('id'),{{'tau',1}},{double(seckey)},{'tau','tau~'},eye(obj.sector_multiplets(seckey)));
                end
                obj.fermionicities('ph') = 1;
                obj.op_is_scalar('ph') = true;
                obj.eta_legs('ph') = {};
            end
end

