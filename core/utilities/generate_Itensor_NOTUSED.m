function Itensor = generate_Itensor(sym,progress_info)
% GENERATE_ITENSOR Generates Dtensor for Symmetry or SuperSym classes.
% If progress_info == true >> a nice progress bar is shown.
            if nargin == 1
                progress_info = false;
            elseif nargin > 2
                error('nargin must be 1 or 2') 
            end
            

            Itensor = 0;
            
            if isa(sym,'Symmetry') && ~isa(sym,'SuperSym')
                if progress_info
                    f = waitbar(0.0,'Itensor calculation progress');
                    totrep = length(sym.CGtensor.get_irrep_set({'m1',1}));
                    i = 1;
                end
                for looprep = sym.CGtensor.get_irrep_set({'m1',1})
                    CG1trunc = sym.CGtensor.irrep_trunc_copy({{{'m1',1},{looprep{1}}}});
                    tmp = dot(CG1trunc,conj(CG1trunc),{'m1'},{'m1'},...
                              {{{'m2','mu'},{'M','mright'}},...
                               {{'m2','M'},{'M','mleft~'},{'alpha','omega_in'}}});
                    
                    tmp = dot(tmp,conj(CG1trunc),{'mleft~','mu'},{'m1','m2'},...   %WROOOOOOOOOOOOOOOOOOOOONG
                              {{},{{'M','mright~'},{'alpha','alpha~'}}});
                    tmp = dot(tmp,sym.one_per_dim,{'mright~'},{'M'},{{},{{'M~','mright~'}}});
                    Itensor = Itensor + dot(tmp,sym.CGtensor,{'mright','M','mright~'},{'m1','m2','M'},...
                                  {{},{{'alpha','omega_out'}}});
                    if progress_info
                        waitbar(i/totrep,f,['Itensor calculation progress: ', num2str(i), '/', num2str(totrep)]);
                        i = i + 1;
                    end
                end
            elseif isa(sym,'SuperSym')
                if progress_info
                    f = waitbar(0.0,'Itensor calculation progress (SUPERSYM)');
                    totrep = 0;
                    for symID = 1:sym.CGtensor.no_of_symmetries
                        totrep = totrep + length(sym.CGtensor.sub_tensors{symID}.get_irrep_set({'m1',1}));
                    end
                    i = 1;
                end
                Itensors = cell(1,sym.CGtensor.no_of_symmetries);
                for symID = 1:sym.CGtensor.no_of_symmetries
                    Itensors{symID} = 0;

                    for looprep = sym.CGtensor.sub_tensors{symID}.get_irrep_set({'m1',1})
                        CG1trunc = sym.CGtensor.sub_tensors{symID}.irrep_trunc_copy({{{'m1',1},{looprep{1}}}});
                        tmp = dot(CG1trunc,conj(CG1trunc),{'m1'},{'m1'},...
                                  {{{'m2','mu'},{'M','mright'}},...
                                   {{'m2','M'},{'M','mleft~'},{'alpha','omega_in'}}});

                        tmp = dot(tmp,conj(sym.CGtensor.sub_tensors{symID}),{'mleft~','mu'},{'m1','m2'},...
                                  {{},{{'M','mright~'},{'alpha','alpha~'}}});
                        tmp = dot(tmp,sym.one_per_dim.sub_tensors{symID},{'mright~'},{'M'},{{},{{'M~','mright~'}}});
                        Itensors{symID} = Itensors{symID} + dot(tmp,sym.CGtensor.sub_tensors{symID},{'mright','M','mright~'},{'m1','m2','M'},...
                                      {{},{{'alpha','omega_out'}}});
                        if progress_info
                            waitbar(i/totrep,f,['Itensor calculation progress (SUPERSYM): ', num2str(i), '/', num2str(totrep)]);
                            i = i + 1;
                        end
                    end
                    if symID == 1
                    elseif symID == 2
                        Itensor = irrep_kron(Itensors{1},Itensors{2},'LNAtensor');
                    else
                        Itensor = irrep_kron(Itensor,Itensors{symID}); 
                    end
                end
            else
                error('sym must be Symmetry or SuperSym object.')
            end
            close(f)
end


