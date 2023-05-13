function out = LSMPS_expval_local_scalar_operator(obj, pos, scalar_op)
%LSMPS_EXPVAL_LOCAL_OPERATOR returns the expectation value of a local
%scalar operator. Works only for scalar operator (with 'tau' and 'tau~' legs only)
    tmp = NTdot(obj.matrices{pos}, obj.schmidt_list{pos}, {'t_right'}, {'t_left'});
    tmp2 = NTdot(tmp,scalar_op,{'tau'},{'tau'});
    out = NTdot(tmp2,NTconj(tmp),{'t_left','tau~','alpha','t_right'},...
                {'t_left','tau','alpha','t_right'});
end

