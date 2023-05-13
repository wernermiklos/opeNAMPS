function envelope = SYM_generate_TEBD_envelope(obj,sitereps1, bondreps, sitereps2)
            % Generates TEBD envelope (4 Clebsch sandwich for the two-site
            % evolver).
            if nargin == 1
                sitereps1 = NTget_irrep_set(obj.CGtensor, {'m2',1});
                sitereps2 = sitereps1 ;
                bondreps = NTget_irrep_set(obj.CGtensor, {'m1',1});
            elseif nargin == 2
                sitereps2 = sitereps1;
                bondreps = NTget_irrep_set(obj.CGtensor, {'m1',1});
            elseif nargin == 3
                sitereps2 = sitereps1;    %otherwise they are different
            elseif nargin > 4
                error('nargin must be 1,2, or 3.')
            end
            CG1 = NTirrep_trunc_copy(obj.CGtensor, {{{'m1',1},bondreps},{{'m2',1},sitereps1}});
            CG2 = NTirrep_trunc_copy(obj.CGtensor, {{{'m1',1},NTget_irrep_set(CG1,{'M',1})},{{'m2',1},sitereps2}});
            envelope = 0;
            for rep = bondreps
                CGtrunc = NTirrep_trunc_copy(CG1, {{{'m1',1},rep}});
                tmp = NTdot(CGtrunc,CG2,{'M'},{'m1'},...
                          {{{'m1','m_left'},{'m2','mu1'},{'alpha','alpha1'}},...
                           {{'m2','mu2'},{'alpha','alpha2'},{'M','m_right'}}});
                tmp2 = NTdot(NTconj(tmp),obj.one_per_dim,{'m_right'},{'M'},{{},{{'M~','m_right'}}});
                tmp2 = NTprime_all_names(tmp2);
                envelope = NTadd(envelope, NTdot(tmp,tmp2,{'m_left','m_right'},{'m_left~','m_right~'}));
            end

end

