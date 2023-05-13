function obj = SITE_dot_define(obj_in, name1, name2, newname,symmetries,kept_irreps)
            % Defines a new operator from operators name1 and name2.
            % result is op(name1) * op(name2) (* |ket>)
            if nargin == 5
                kept_irreps = {};
            elseif nargin > 6
                error('dot_define has either 5 or 6 arguments.')
            end
            obj = obj_in;
            obj.fermionicities(newname) = obj.fermionicities(name1)*obj.fermionicities(name2);
            etanum1 = 0;
            etanum2 = 0;
            for legname1 = obj.operators(name1).leg_names
                if strcmp(legname1{1}(1:3),'eta')
                    etanum1 = etanum1 + 1; 
                end
            end
            for legname2 = obj.operators(name2).leg_names
                if strcmp(legname2{1}(1:3),'eta')
                   etanum2 = etanum2 + 1; 
                end
            end
            if etanum2 > 0
                error('op2 cannot have eta legs.')
            end
            if obj.op_is_scalar(name1)
                if obj.op_is_scalar(name2)
                    obj.op_is_scalar(newname) = true;
                end
                obj.op_is_scalar(newname) = false;
                obj.operators(newname) = NTdot(obj.operators(name1),obj.operators(name2),{'tau'},{'tau~'});
                obj.op_specifier_legs(newname) = obj.op_specifier_legs(name1);
            else
                obj.op_is_scalar(newname) = false;
                if obj.op_is_scalar(name2)
                    obj.operators(newname) = NTdot(obj.operators(name1),obj.operators(name2),{'tau'},{'tau~'});
                    obj.op_specifier_legs(newname) = obj.op_specifier_legs(name1);
                else
                    oprepset1 = NTget_irrep_set_1sym(obj.operators(name1), {'omega',2});
                    oprepset2 = NTget_irrep_set_1sym(obj.operators(name2), {'omega',2});
                    siterepset = NTget_irrep_set_1sym(obj.operators('id'),{'tau', 1});
                    for symID = 1:length(symmetries)
                        CG1 = NTirrep_trunc_copy(symmetries{symID}.CGtensor, {{{'M',1},siterepset{symID}},{{'m2',1},oprepset1{symID}},{{'m1',1},siterepset{symID}}});
                        CG2 = NTirrep_trunc_copy(symmetries{symID}.CGtensor, {{{'M',1},siterepset{symID}},{{'m2',1},oprepset2{symID}},{{'m1',1},siterepset{symID}}});
                        CG3 = NTirrep_trunc_copy(symmetries{symID}.CGtensor, {{{'m2',1},oprepset2{symID}},{{'m1',1},oprepset1{symID}}});
                        CG4 = NTirrep_trunc_copy(symmetries{symID}.CGtensor, {{{'M',1},siterepset{symID}},{{'m1',1},siterepset{symID}}});
                        CG4 = NTdot(CG4,symmetries{symID}.one_per_dim,{'M'},{'M~'});
                        tmp = NTdot(NTconj(CG1),NTconj(CG2),{'m1'},{'M'},...
                                  {{{'alpha','omega1'},{'M','mu~'},{'m2','M_op1'}},...
                                   {{'alpha','omega2'},{'m1','mu'},{'m2','M_op2'}}});
                        tmp = NTdot(tmp,CG3,{'M_op1','M_op2'},{'m1','m2'},...
                                  {{},{{'alpha','eta'},{'M','M_op'}}});
                        tmp = NTdot(tmp,CG4,{'mu','M_op','mu~'},{'m1','m2','M'},{{},{{'alpha','omega'}}});
                        if symID == 1
                            con_ten = NTcopy(tmp);
                        else
                            con_ten = NTirrep_kron(con_ten,tmp,'LNAtensor');
                        end
                    end
                    tmp = NTdot(obj.operators(name1),obj.operators(name2),{'tau'},{'tau~'},{{{'omega','omega1'}},{{'omega','omega2'}}});
                    newetalegname = ['eta',num2str(etanum1+etanum2+1)];
                    obj.operators(newname) = NTdot(tmp,con_ten,{'omega1','omega2'},{'omega1','omega2'},...
                                                 {{},{{'eta',newetalegname}}});
                    if ~isempty(kept_irreps)
                        obj.operators(newname) = NTirrep_trunc_copy(obj.operators(newname),{{{'omega',2},kept_irreps}});
                        if length(kept_irreps) == 1 && all(kept_irreps{1} == 1)    % Operator is scalar. We remove the omega leg.
                            omega_sectors = NTget_leg_sectors(obj.operators(newname),'omega');
                            omega_remover = NAtensor({'omega'},{'i'},{[1,2,3]},obj.no_of_symmetries);
                            for secID = 1:length(omega_sectors)
                                 omega_remover = NTset_block(omega_remover,{{'omega',1},{'omega',2},{'omega',3}},omega_sectors{secID}{1},{'omega'},[1]);
                            end
                            obj.operators(newname) = NTdot(obj.operators(newname),omega_remover,{'omega'},{'omega'});
                            obj.operators(newname) = NTmatch_irreps(obj.operators(newname),{'tau',1},{'tau~',1});
                            obj.op_is_scalar(newname) = true;
                        end
                    end
                    obj.operators(newname) = NTkill_small_blocks(obj.operators(newname), 1.e-14);
                    obj.eta_legs(newname) = [obj.eta_legs(name1),newetalegname];
                end
            end
           

end

