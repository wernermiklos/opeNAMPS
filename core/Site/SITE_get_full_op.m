function out = SITE_get_full_op(obj,name,symmetries)
            if obj.op_is_scalar(name)
                delta_mu = NAtensor({'mu','mu~'},{'o','i'},{[1],[1]},obj.no_of_symmetries);
                for seckey = obj.sector_multiplets.keys
                    delta_mu = NTset_block(delta_mu,{{'mu',1}},{double(seckey{1})},{'mu','mu~'},eye(obj.multiplet_dims(seckey{1})));
                end
                out = NTsimple_mult(obj.operators(name),delta_mu,{{{'tau',1},{'mu',1}}});
            else
                siteirreps = SITE_get_site_irreps(obj);
                if length(symmetries) == 1
                    CG = symmetries{1}.CGtensor;
                    CG = NTirrep_trunc_copy(CG, {{{'m1',1},siteirreps},{{'M',1},siteirreps}});
                    out = NTdot(obj.operators(name),NTconj(CG),{'omega'},{'alpha'},{{},{{'m1','mu'},{'M','mu~'},{'m2','mu_op'}}});
                else
                    CG = NTirrep_kron(symmetries{1}.CGtensor, symmetries{2}.CGtensor, 'LNAtensor');
                    for i = 3:length(symmetries)
                        CG = LNTirrep_kron(CG, symmetries{i}.CGtensor); 
                    end
                    CG = LNTirrep_trunc_copy(CG, {{{'m1',1},siteirreps},{{'M',1},siteirreps}});
                    out = LNTdot(obj.operators(name),LNTconj(CG),{'omega'},{'alpha'},{{},{{'m1','mu'},{'M','mu~'},{'m2','mu_op'}}});
                end 
            end  

end

