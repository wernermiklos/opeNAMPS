function obj = SITE_define_tensor_operator(obj, ...
                                        name, ...                  % char array, the name of the operator (like 'fdag')
                                        symmetries, ...            % cell array of symmetries
                                        opirrep, ...               % The representation ("order") to which the operator multiplet belongs
                                        oplist, ...                % List of operators (the list index corresponds to the "M_op" index of the tensor operator)
                                        varargin)                 % Tolerance value (10^(-12) by default). It is used to check if the op is indeed a tensor op.
             parameters = struct();
             % Default parameters
             parameters.FERMIONICITY = +1;
             parameters.TOLERANCE = 1.e-12;
                % Update according to varargin
             parameters = parameter_updater(parameters,varargin);
             
             fermionicity = parameters.FERMIONICITY;
             tolerance = parameters.TOLERANCE;
             
             super_one_per_dim = symmetries{1}.one_per_dim;
             super_CGtensor = symmetries{1}.CGtensor;
             if obj.no_of_symmetries == 1
                 CGtrunc = NTirrep_trunc_copy(super_CGtensor,{{{'m1',1},cellfun(@(x) double(x), obj.sector_multiplets.keys, 'UniformOutput', false)}, ...
                                                        {{'m2',1},{opirrep}}, ...
                                                        {{'M',1},cellfun(@(x) double(x), obj.sector_multiplets.keys, 'UniformOutput', false)}});
             else
                 for symID = 2:(obj.no_of_symmetries)
                     if symID == 2
                         super_one_per_dim = NTirrep_kron(super_one_per_dim, symmetries{symID}.one_per_dim, 'LNAtensor');
                         super_CGtensor = NTirrep_kron(super_CGtensor, symmetries{symID}.CGtensor, 'LNAtensor');
                     else
                         super_one_per_dim = LNTirrep_kron(super_one_per_dim, symmetries{symID}.one_per_dim);
                         super_CGtensor = LNTirrep_kron(super_CGtensor, symmetries{symID}.CGtensor);
                     end
                 end
             CGtrunc = LNTirrep_trunc_copy(super_CGtensor,{{{'m1',1},cellfun(@(x) double(x), obj.sector_multiplets.keys, 'UniformOutput', false)}, ...
                                                        {{'m2',1},{opirrep}}, ...
                                                        {{'M',1},cellfun(@(x) double(x), obj.sector_multiplets.keys, 'UniformOutput', false)}});
             end

             
             obj.fermionicities(name) = fermionicity;

             if all(opirrep == 1)   % If it's a scalar operator
                 if length(oplist) > 1
                     error('Scalar operators are one-dimensional.')
                 end
                 obj.op_is_scalar(name) = true;
                 obj.operators(name) = NTdot(oplist{1},super_one_per_dim,{'mu','mu~'},{'M~','M'});
                 op_copy = NTcopy(oplist{1});
                 op_copy = NTmatch_irreps(op_copy,{'tau',1},{'tau~',1});
                 delta_mu = NAtensor({'mu','mu~'},{'o','i'},{[1],[1]},obj.no_of_symmetries);
                 for seckey = obj.sector_multiplets.keys
                     delta_mu = NTset_block(delta_mu,{{'mu',1}},{double(seckey{1})},{'mu','mu~'},eye(obj.multiplet_dims(seckey{1})));
                 end
                 test = NTsimple_mult(obj.operators(name),delta_mu,{{{'tau',1},{'mu',1}}});
                 test = NTsubtr(test,op_copy);

                 if NTget_max_tensor_element(test) > tolerance
                     warning(['Difference between original and stored op: ', num2str(NTget_max_tensor_element(test))]);
                 end
             else
                 irrepdim = 1;
                for symID = 1:obj.no_of_symmetries
                    irrepdim = irrepdim * symmetries{symID}.irrep_dimensions(opirrep(symID));  
                end
                 if length(oplist) ~= irrepdim
                     error(['Operator of irrep [', num2str(opirrep), '] should be ', num2str(irrepdim), ' dimensional.'])
                 end
                 obj.op_is_scalar(name) = false;
                 tensorop = NAtensor({'tau','mu','tau~','mu~','M_op'},{'o','o','i','i','o'},{[1],[1],[2],[2],[3]},obj.no_of_symmetries);
                 [irreps,blocks,shapes] = NTexport_all_data_blocks(oplist{1},{{'tau',1},{'tau~',1}},{'tau','mu','tau~','mu~'});
                 for bID = 1:length(irreps)
                     newirreps = [irreps{bID},{opirrep}];
                     newblock = zeros(prod(shapes{bID}),irrepdim);
                     newblock(:,1) = blocks{bID}(:);
                     for M_op = 2:irrepdim
                         block_tmp = NTget_block(oplist{M_op},{{'tau',1},{'tau~',1}},irreps{bID},{'tau','mu','tau~','mu~'});
                         newblock(:,M_op) = block_tmp(:);
                     end
                     newblock = reshape(newblock,[shapes{bID},irrepdim]);
                     tensorop = NTset_block(tensorop, {{'tau',1},{'tau~',1},{'M_op',1}},newirreps,{'tau','mu','tau~','mu~','M_op'},newblock);
                 end
                 if obj.no_of_symmetries == 1
                    CGtmp = NTdot(CGtrunc,super_one_per_dim,{'M'},{'M~'});
                 
                    obj.operators(name) = NTdot(tensorop,CGtmp,{'mu','mu~','M_op'},{'m1','M','m2'},{{},{{'alpha','omega'}}});
                 
                    test = NTdot(obj.operators(name), NTconj(CGtrunc),{'omega'},{'alpha'},{{},{{'m1','mu'},{'M','mu~'},{'m2','M_op'}}});
                 else
                    CGtmp = LNTdot(CGtrunc,super_one_per_dim,{'M'},{'M~'});
                 
                    obj.operators(name) = LNTdot(tensorop,CGtmp,{'mu','mu~','M_op'},{'m1','M','m2'},{{},{{'alpha','omega'}}});
                 
                    test = LNTdot(obj.operators(name), LNTconj(CGtrunc),{'omega'},{'alpha'},{{},{{'m1','mu'},{'M','mu~'},{'m2','M_op'}}});
                 end
                 
                 test = NTsubtr(test,tensorop);

                 if NTget_max_tensor_element(test) > tolerance
                     warning(['Difference between original and stored op: ', num2str(NTget_max_tensor_element(test))]);
                 end
             end
             obj.eta_legs(name) = {};
end

