function out = SYM_generate_finitesegmentDMdiagonalizer(obj,seglength,sitereps,bondreps,bondreps_left)
%SYM_GENERATE_FINITESEGMENTDMDIAGONALIZER Summary of this function goes here
%   Detailed explanation goes here
    if nargin == 2
                sitereps = NTget_irrep_set(obj.CGtensor,{'m2',1});
                bondreps = NTget_irrep_set(obj.CGtensor,{'m1',1});
                bondreps_left = bondreps;
    elseif nargin == 3
                bondreps = NTget_irrep_set(obj.CGtensor,{'m1',1});
                bondreps_left = bondreps;
    elseif nargin == 4
                bondreps_left = bondreps;
    end
    CG = NTirrep_trunc_copy(obj.CGtensor,{{{'m1',1},bondreps},{{'m2',1},sitereps},{{'M',1},bondreps}});
    CG2 = NTirrep_trunc_copy(obj.CGtensor,{{{'m2',1},sitereps}});
    sqrt_one_per_dim = NTcopy(obj.one_per_dim);
    for key = sqrt_one_per_dim.data.keys
        sqrt_one_per_dim.data(key{1}) = sqrt(sqrt_one_per_dim.data(key{1}));
    end
    out = 0;
    for repleft = bondreps_left
        disp(repleft)
        tmp = NTrename_legs(NTirrep_trunc_copy(CG,{{{'m1',1},{repleft{1}}}}),...
                       {'m1','m2','M','alpha'},{'m_left','M','m_right','alpha1'});
        for i = 2:seglength
            disp(i);
            if i ~= seglength
                tmp = NTdot(tmp,NTconj(CG2),{'M'},{'m1'},{{},{{'m2','mu'},{'alpha',['beta',num2str(i-1)]}}});
                tmp = NTdot(tmp,CG,{'m_right','mu'},{'m1','m2'},{{},{{'M','m_right'},{'alpha',['alpha',num2str(i)]}}});
            else
                tmp = NTdot(tmp,NTconj(NTdot(CG2,sqrt_one_per_dim,{'M'},{'M~'})), ...
                        {'M'},{'m1'},{{},{{'m2','mu'},{'alpha',['beta',num2str(i-1)]}}});
                tmp = NTdot(tmp,NTdot(CG,sqrt_one_per_dim,{'M'},{'M~'}),...
                        {'m_right','mu'},{'m1','m2'},{{},{{'M','m_right'},{'alpha',['alpha',num2str(i)]}}});
            end
        end
        tmptild = NTprime_all_names(NTconj(tmp));
        out = NTadd(out, NTdot(tmp,tmptild,{'m_left','m_right','M'},{'m_left~','m_right~','M~'}));
    end
end

