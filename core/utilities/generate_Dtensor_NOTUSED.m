function Dtensor = generate_Dtensor(sym,progress_info)
% GENERATE_DTENSOR Generates Dtensor for Symmetry or SuperSym classes.
% If progress_info == true >> a nice progress bar is shown.
            if nargin == 1
                progress_info = false;
            elseif nargin > 2
                error('nargin must be 1 or 2') 
            end
            Dtensor = 0;
            
            if isa(sym,'Symmetry') && ~isa(sym,'SuperSym')
                if progress_info
                    f = waitbar(0.0,'Dtensor calculation progress');
                    totrep = length(sym.CGtensor.get_irrep_set({'m1',1}));
                    i = 1;
                 end
                CG2trunc = sym.CGtensor.irrep_trunc_copy({{{'m1',1},sym.CGtensor.get_irrep_set({'m2',1})}, ...
                                                          {{'M',1}, sym.CGtensor.get_irrep_set({'m2',1})}});

                for looprep = sym.CGtensor.get_irrep_set({'m1',1})
                    CG1trunc = sym.CGtensor.irrep_trunc_copy({{{'m1',1},{looprep{1}}}});
                    tmp = dot(CG1trunc,conj(CG2trunc),{'m2'},{'m1'},...
                              {{{'m1','mleft'},{'M','mright'}},...
                               {{'m2','M'},{'M','mu~'},{'alpha','omega_loc'}}});
                    tmp = dot(tmp,conj(CG1trunc),{'mleft','mu~'},{'m1','m2'},...
                              {{},{{'M','mright~'},{'alpha','alpha~'}}});
                    tmp = dot(tmp,sym.one_per_dim,{'mright~'},{'M'},{{},{{'M~','mright~'}}});
                    Dtensor = Dtensor + dot(tmp,sym.CGtensor,{'mright','M','mright~'},{'m1','m2','M'},...
                                  {{},{{'alpha','omega_out'}}});
                    if progress_info
                        waitbar(i/totrep,f,['Dtensor calculation progress: ', num2str(i), '/', num2str(totrep)]);
                        i = i + 1;
                    end
                end
            elseif isa(sym,'SuperSym')
                if progress_info
                    f = waitbar(0.0,'Dtensor calculation progress (SUPERSYM)');
                    totrep = 0;
                    for symID = 1:sym.CGtensor.no_of_symmetries
                        totrep = totrep + length(sym.CGtensor.sub_tensors{symID}.get_irrep_set({'m1',1}));
                    end
                    i = 1;
                end
                Dtensors = cell(1,sym.CGtensor.no_of_symmetries);
                for symID = 1:sym.CGtensor.no_of_symmetries
                    Dtensors{symID} = 0;
                    CG2trunc = sym.CGtensor.sub_tensors{symID}.irrep_trunc_copy({{{'m1',1},sym.CGtensor.get_irrep_set({'m2',1})}, ...
                                                                                {{'M',1}, sym.CGtensor.get_irrep_set({'m2',1})}});

                    for looprep = sym.CGtensor.sub_tensors{symID}.get_irrep_set({'m1',1})
                        CG1trunc = sym.CGtensor.sub_tensors{symID}.irrep_trunc_copy({{{'m1',1},{looprep{1}}}});
                        tmp = dot(CG1trunc,conj(CG2trunc),{'m2'},{'m1'},...
                                  {{{'m1','mleft'},{'M','mright'}},...
                                   {{'m2','M'},{'M','mu~'},{'alpha','omega_loc'}}});
                        tmp = dot(tmp,conj(CG1trunc),{'mleft','mu~'},{'m1','m2'},...
                                  {{},{{'M','mright~'},{'alpha','alpha~'}}});
                        tmp = dot(tmp,sym.one_per_dim.sub_tensors{symID},{'mright~'},{'M'},{{},{{'M~','mright~'}}});
                        Dtensors{symID} = Dtensors{symID} + dot(tmp,sym.CGtensor.sub_tensors{symID},...
                                                                {'mright','M','mright~'},{'m1','m2','M'},...
                                                                {{},{{'alpha','omega_out'}}});
                        if progress_info
                            waitbar(i/totrep,f,['Dtensor calculation progress (SUPERSYM): ', num2str(i), '/', num2str(totrep)]);
                            i = i + 1;
                        end
                    end
                    if symID == 1
                    elseif symID == 2
                        Dtensor = irrep_kron(Dtensors{1},Dtensors{2},'LNAtensor');
                    else
                        Dtensor = irrep_kron(Dtensor,Dtensors{symID}); 
                    end
                end
            else
                error('sym must be Symmetry or SuperSym object.')
            end
            close(f)
end



