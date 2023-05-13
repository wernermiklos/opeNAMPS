function out = ABMPS_measure_singlesite_op(obj, oppos, op)
%ABMPS_MEASURE_SINGLESITEOP Summary of this function goes here
%   Detailed explanation goes here
    if ~ABMPS_is_left(obj)
        error('MPS should be left canonical!!')
    end
    if oppos == 1
        lefttmp = NTdot(obj.left_matrices{1},op,{'tau'},{'tau~'});
        lefttmp = NTdot(lefttmp,NTconj(obj.left_matrices{1}),{'tau~','t_in'},{'tau','t_in'},{{},{{'t_out','t_out~'}}});
        for pos = 2:(obj.chain_length)
           lefttmp = NTdot(lefttmp,obj.left_matrices{pos},{'t_out'},{'t_in'});
           lefttmp = NTdot(lefttmp,NTconj(obj.left_matrices{pos}),{'t_out~','tau'},{'t_in','tau'},{{},{{'t_out','t_out~'}}});
        end
        lefttmp = NTdot(lefttmp,obj.schmidt,{'t_out'},{'t_left'});
        out = NTdot(lefttmp,NTconj(obj.schmidt),{'t_right','t_out~'},{'t_right','t_left'});
    else
        lefttmp = NTdot(obj.left_matrices{1},NTconj(obj.left_matrices{1}),{'tau','t_in'},{'tau','t_in'},{{},{{'t_out','t_out~'}}});
        for pos = 2:(obj.chain_length)
           lefttmp = NTdot(lefttmp,obj.left_matrices{pos},{'t_out'},{'t_in'});
           if pos == oppos
                lefttmp = NTdot(lefttmp,op,{'tau'},{'tau'},{{},{{'tau~','tau'}}});
           end
           lefttmp = NTdot(lefttmp,NTconj(obj.left_matrices{pos}),{'t_out~','tau'},{'t_in','tau'},{{},{{'t_out','t_out~'}}});
        end
        lefttmp = NTdot(lefttmp,obj.schmidt,{'t_out'},{'t_left'});
        out = NTdot(lefttmp,NTconj(obj.schmidt),{'t_right','t_out~'},{'t_right','t_left'});
    end
end

