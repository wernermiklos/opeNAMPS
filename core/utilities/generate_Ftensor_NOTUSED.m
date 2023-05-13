function Ftensor = generate_Ftensor(sym, progress_info)
%GENERATE_FTENSOR Summary of this function goes here
%   Detailed explanation goes here
    if nargin == 1
        progress_info = false;
    elseif nargin > 2
        error('nargin must be 1 or 2') 
    end
    Ftensor = 0;
            
    if isa(sym,'Symmetry') && ~isa(sym,'SuperSym')
        if progress_info
            f = waitbar(0.0,'Ftensor calculation progress');
            totrep = length(sym.CGtensor.get_irrep_set({'m1',1}))*length(sym.CGtensor.get_irrep_set({'m2',1}));
            i = 1;
        end
            CG2trunc = sym.CGtensor.irrep_trunc_copy({{{'m1',1},sym.CGtensor.get_irrep_set({'m2',1})}, ...
                                                      {{'M',1}, sym.CGtensor.get_irrep_set({'m2',1})}});
                for looprep = sym.CGtensor.get_irrep_set({'m1',1})
                    for oprep = sym.CGtensor.get_irrep_set({'m2',1})
                        CG1trunc = sym.CGtensor.irrep_trunc_copy({{{'m1',1},{looprep{1}}}});
                        CG3trunc = CG1trunc.irrep_trunc_copy({{{'m2',1},{oprep{1}}}});
                        CG4trunc = CG2trunc.irrep_trunc_copy({{{'m1',1},{oprep{1}}}});
                        CG5trunc = sym.CGtensor.irrep_trunc_copy({{{'m1',1},CG1trunc.get_irrep_set({'M',1})},...
                                                                  {{'m2',1},CG4trunc.get_irrep_set({'M',1})}});
                        CG6trunc = sym.CGtensor.irrep_trunc_copy({{{'m1',1},CG3trunc.get_irrep_set({'M',1})}});                                      
                        tmp1 = dot(CG5trunc,sym.one_per_dim,...
                                                {'M'},{'M~'}, ...
                                                {{{'m1','mright'},{'m2','Mout'},{'alpha','omega_out'}},...
                                                 {{'M','mright~'}}});
                        tmp1 = dot(tmp1, CG1trunc,...
                                        {'mright'},{'M'},...
                                        {{},{{'m1','mleft'},{'m2','mu'}}});
                        tmp2 = dot(conj(CG3trunc),conj(CG6trunc),...
                                   {'M'},{'m1'},...
                                   {{{'alpha','omega_in'},{'m1','mleft'},{'m2','Min'}},...
                                    {{'alpha','alpha~'},{'m2','mu~'},{'M','mright~'}}});
                        tmp12 = dot(tmp1,tmp2,{'mleft','mright~'},{'mleft','mright~'});      
                        tmp3 = dot(conj(CG2trunc),CG4trunc,...
                                       {'m2'},{'m2'},...
                                       {{{'m1','mu'},{'alpha','omega_loc'},{'M','mu~'}},{{'m1','Min'},{'M','Mout'},{'alpha','xi'}}});
                        Ftensor = Ftensor + dot(tmp12,tmp3,{'mu','mu~','Min', 'Mout'},{'mu','mu~','Min','Mout'});
                        if progress_info
                            waitbar(i/totrep,f,['Ftensor calculation progress: ', num2str(i), '/', num2str(totrep)]);
                            i = i + 1;
                        end
                    end
                end
            elseif isa(sym,'SuperSym')
                if progress_info
                    f = waitbar(0.0,'Dtensor calculation progress (SUPERSYM)');
                    totrep = 0;
                    for symID = 1:sym.CGtensor.no_of_symmetries
                        totrep = totrep + length(sym.CGtensor.sub_tensors{symID}.get_irrep_set({'m1',1}))*length(sym.CGtensor.sub_tensors{symID}.get_irrep_set({'m2',1}));
                    end
                    i = 1;
                end
                Ftensors = cell(1,sym.CGtensor.no_of_symmetries);
                for symID = 1:sym.CGtensor.no_of_symmetries
                    Ftensors{symID} = 0;
                    CG2trunc = sym.CGtensor.sub_tensors{symID}.irrep_trunc_copy({{{'m1',1},sym.CGtensor.get_irrep_set({'m2',1})}, ...
                                                                                {{'M',1}, sym.CGtensor.get_irrep_set({'m2',1})}});

                    for looprep = sym.CGtensor.sub_tensors{symID}.get_irrep_set({'m1',1})
                        CG1trunc = sym.CGtensor.sub_tensors{symID}.irrep_trunc_copy({{{'m1',1},{looprep{1}}}});
                        tmp1 = dot(sym.CGtensor.sub_tensors{symID},sym.one_per_dim.sub_tensors{symID},...
                                            {'M'},{'M~'}, ...
                                            {{{'m1','mright'},{'m2','Mout'},{'alpha','omega_out'}},...
                                             {{'M','mright~'}}});
                        tmp1 = dot(tmp1, CG1trunc,...
                                        {'mright'},{'M'},...
                                        {{},{{'m1','mleft'},{'m2','mu'}}});
                        tmp2 = dot(conj(CG1trunc),conj(sym.CGtensor.sub_tensors{symID}),...
                                   {'M'},{'m1'},...
                                   {{{'alpha','omega_in'},{'m1','mleft'},{'m2','Min'}},...
                                    {{'alpha','alpha~'},{'m2','mu~'},{'M','mright~'}}});
                        tmp12 = dot(tmp1,tmp2,{'mleft','mright~'},{'mleft','mright~'});      
                        tmp3 = dot(conj(CG2trunc),CG2trunc,...
                                   {'m2'},{'m2'},...
                                   {{{'m1','mu'},{'alpha','omega_loc'},{'M','mu~'}},{{'m1','Min'},{'M','Mout'},{'alpha','xi'}}});
                        Ftensors{symID} = Ftensors{symID} + dot(tmp12,tmp3,{'mu','mu~','Min', 'Mout'},{'mu','mu~','Min','Mout'});
                        if progress_info
                            waitbar(i/totrep,f,['Dtensor calculation progress (SUPERSYM): ', num2str(i), '/', num2str(totrep)]);
                            i = i + 1;
                        end
                    end
                    if symID == 1
                    elseif symID == 2
                        Ftensor = irrep_kron(Ftensors{1},Ftensors{2},'LNAtensor');
                    else
                        Ftensor = irrep_kron(Ftensor,Ftensors{symID}); 
                    end
                end
            else
                error('sym must be Symmetry or SuperSym object.')
            end
            close(f)
end

