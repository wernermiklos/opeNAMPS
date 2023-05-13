function obj = SYM_create(CGtensor,info)
            obj.CGtensor = NAtensor({'m1','m2','M','alpha'}, ...  
                                    {'i','i','o','o'}, ...
                                    {[1],[2],[3],[1,2,3]}, ...
                                    1);                               % We create an empty tensor 
            obj.CGtensor = NTadd(obj.CGtensor,CGtensor);                   % We fill it with values   !!!!!
            if nargin == 1
                info = '';
            elseif nargin > 2
                error('nargin must be 1 or 2');
            end
            if ~isa(info,'char')
                error('info must be char (i.e. characters between apostrophes)')
            end
            obj.type = 'Symmetry';
            obj.info = info;
            obj.irrep_dimensions = containers.Map('KeyType','double','ValueType','double');
            obj.active_reps = {[],[],[]};
            for key = obj.CGtensor.data.keys()
                keyvalues = double(key{1});
                datablock = obj.CGtensor.data(key{1});
                shape = [size(datablock), ones(1,4-ndims(datablock))];
                for i = 1:length(keyvalues)
                    obj.active_reps{i}(end+1) = keyvalues(i);
                    if obj.irrep_dimensions.isKey(keyvalues(i))
                        if shape(i) ~= obj.irrep_dimensions(keyvalues(i))     % Here we exploit the known leg ordering of the CGtensor
                            error('CGtensor is corrupted. Irrep dimensions are inconsistent.')
                        end
                    else
                        obj.irrep_dimensions(keyvalues(i)) = shape(i);
                    end
                end
            end
            for i = 1:3
                obj.active_reps{i} = sort(unique(obj.active_reps{i}));
            end
            % Some basic functions/storages (these can be changed for specific
            % symmetries in order to improve performance)
            obj.selection_rule = @(sym,Gamma1, Gamma2, Gamma) SYM_selection_rule_basic(sym,Gamma1,Gamma2,Gamma);
            obj.fusion_rule = @(sym,Gamma1,Gamma2) SYM_fusion_rule_basic(sym,Gamma1,Gamma2);
            obj.conj_reps = containers.Map('KeyType','double','ValueType','double');
            for Gamma1 = obj.active_reps{1}
                CGtr = NTirrep_trunc_copy(obj.CGtensor,{{{'m1',1},{Gamma1}},{{'M',1},{[1]}}});
                secs = NTget_leg_sectors(CGtr,'m2');
                if ~isempty(secs)
                    obj.conj_reps(Gamma1) = secs{1}{1}{1};
                end
            end

            
            obj.one_per_dim = NAtensor({'M~','M'},...
                                       {'i','o'},...
                                       {[1],[1]},...
                                       1);
            for key = obj.irrep_dimensions.keys
                obj.one_per_dim = NTset_block(obj.one_per_dim, {{'M',1}},{key{1}},{'M~','M'},...
                                          1.0/obj.irrep_dimensions(key{1}) * eye( obj.irrep_dimensions(key{1})));
            end

            obj.Qnum_to_Gamma = @(x) x;
            obj.Gamma_to_Qnum = @(x) x;


end

