function PANDA = SYM_generate_PANDAtensor(obj)
%SYM_GENERATE_PANDATENSOR Summary of this function goes here
%   Detailed explanation goes here
PANDA_projector = NAtensor({'M','alpha'},{'i','i'},{[3],[1,2,3]},1);
PANDA_tmp = NTirrep_trunc_copy(obj.CGtensor,{{{'M',1},{1}}});
active_blocks = NTget_active_irrep_values(PANDA_tmp,{{'alpha',1},{'alpha',2},{'alpha',3}});
for bID = 1:length(active_blocks)
    PANDA_projector = NTset_block(PANDA_projector,{{'alpha',1},{'alpha',2},{'alpha',3}},...
        active_blocks{bID},{'M','alpha'},[1]);
end
PANDA = NTdot(PANDA_tmp,PANDA_projector,{'alpha','M'},{'alpha','M'});


end

