function out = LSMPS_expval_scalar_bond_operator(obj, pos, Twosite_red_op)
%LSMPS_EXPVAL_LOCAL_OPERATOR returns the expectation value of a two site 
%bond (nearest neighbor) operator. Works correctly only for scalar
%operator. 
    pos1 = pos;
    pos2 = mod(pos,obj.chain_length) + 1;
    tmp = NTdot(obj.matrices{pos1},obj.matrices{pos2},{'t_right'},{'t_left'},...
                {{{'tau','tau_1'},{'alpha','alpha_1'}}, ...
                 {{'tau','tau_2'},{'alpha','alpha_2'}}});
    psi2site = NTdot(tmp, obj.schmidt_list{pos2}, {'t_right'}, {'t_left'});
    tmp = NTdot(psi2site,Twosite_red_op, ...
                {'tau_1','alpha_1','tau_2','alpha_2'}, ...
                {'tau_1','alpha_1','tau_2','alpha_2'});
    out = NTdot(tmp,NTconj(psi2site),{'t_left','tau_1~','alpha_1~', 'tau_2~', 'alpha_2~', 't_right'},...
                {'t_left','tau_1','alpha_1', 'tau_2', 'alpha_2', 't_right'});
end

