function out = ABMPS_scalarprod(MPS1, MPS2)
%ABMPS_SCALARPROD Calculates the scalar product of MPS1 and MPS2. MPS1 is
%   conjugated. <MPS1|MPS2>
%   Detailed explanation goes here
    if ~ABMPS_is_left(MPS1) ||  ~ABMPS_is_left(MPS2) 
        error('MPS1 & MPS2 should be left canonical!!')
    end
    if MPS1.chain_length ~= MPS2.chain_length
        error('Chain length should be equal for MPS1 and MPS2');
    end
    lefttmp = NTdot(MPS2.left_matrices{1},NTconj(MPS1.left_matrices{1}),{'tau','t_in'},{'tau','t_in'},{{},{{'t_out','t_out~'}}});
    for pos = 2:(MPS1.chain_length)
       lefttmp = NTdot(lefttmp,MPS2.left_matrices{pos},{'t_out'},{'t_in'});
       lefttmp = NTdot(lefttmp,NTconj(MPS1.left_matrices{pos}),{'t_out~','tau'},{'t_in','tau'},{{},{{'t_out','t_out~'}}});
    end
    lefttmp = NTdot(lefttmp,MPS2.schmidt,{'t_out'},{'t_left'});
    out = NTdot(lefttmp,NTconj(MPS1.schmidt),{'t_right','t_out~'},{'t_right','t_left'});

end

