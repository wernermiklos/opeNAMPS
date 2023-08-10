function out = COUP_generate_TwoSite_full_OpProd(symmetries, site1, site2, opname1, opname2, cp_combinations )
%GENERATE_TWOSITE_OPPROD Generates the two site Operator Product for TEBD.
%   op1 is the left site, op2 is the right, and the operator is op1 * op2!
%   (left-to-right ordering is important)
%   cp_combinations format {{[m1,m2],val},{[m1,m2],val}};
    
    
    if ~isempty(site1.eta_legs(opname1)) || ~isempty(site2.eta_legs(opname2))
        error('not implemented for combined operators with eta legs');
    end
    
    if site1.fermionicities(opname1)*site2.fermionicities(opname2) < 0
        error('total fermionicity should be +1')
    end
    
    op1 = SITE_get_full_op(site1, opname1, symmetries);
    op2 = SITE_get_full_op(site2, opname2, symmetries);
    
    if site1.fermionicities(opname1) < 0
       op1 = NTdot(site1.operators('ph'),op1,{'tau~'},{'tau'});
    end

    
 
    
    
    op1 = NTrename_legs(op1, {'mu','mu~','tau','tau~'},{'mu_1','mu_1~','tau_1','tau_1~'});
    op2 = NTrename_legs(op2, {'mu','mu~','tau','tau~'},{'mu_2','mu_2~','tau_2','tau_2~'});
    
    if site1.op_is_scalar(opname1) && site2.op_is_scalar(opname2)
        if length(cp_combinations) > 1 || ~isequal(cp_combinations{1}{1},[1,1])
            error('For two scalar operators cp_combination = {{[1,1],value}}');
        end
        out = NTmult(NTsimple_mult(op1,op2), cp_combinations{1}{2}); 
        out = NTsplit_irrep(out,{'tau_1',1},{'tau_1~','mu_1~'});    % Scalar ops are irrep diagonal, but the generic two-site op is not.
        out = NTsplit_irrep(out,{'tau_2',1},{'tau_2~','mu_2~'});    % Scalar ops are irrep diagonal, but the generic two-site op is not.
        
    elseif ~site1.op_is_scalar(opname1) && ~site2.op_is_scalar(opname2)
        mu_op1_sectors = NTget_leg_sectors(op1, 'mu_op');
        mu_op2_sectors = NTget_leg_sectors(op2, 'mu_op');
        if length(mu_op1_sectors)>1 || length(mu_op2_sectors)>1
            error('Something is wrong. Op1 or op2 is not a simple tensor operator.');
        end
        mu_op1_dim = mu_op1_sectors{1}{2};
        op1_irrep = mu_op1_sectors{1}{1}{1};
        mu_op2_dim = mu_op2_sectors{1}{2};
        op2_irrep = mu_op2_sectors{1}{1}{1};


        cp_tensor = NAtensor({'mu_op1','mu_op2'},{'i','i'},{[1],[2]},length(symmetries));
        cpblock = zeros(mu_op1_dim, mu_op2_dim);
        for i = 1:length(cp_combinations)
            cpblock(cp_combinations{i}{1}(1),cp_combinations{i}{1}(2)) = cp_combinations{i}{2};
        end
        cp_tensor = NTset_block(cp_tensor, {{'mu_op1',1},{'mu_op2',1}},{op1_irrep,op2_irrep},{'mu_op1','mu_op2'},cpblock);
        out = NTdot(op1,cp_tensor,{'mu_op'},{'mu_op1'});
        out = NTdot(out,op2,{'mu_op2'},{'mu_op'});
    else
        error('Scalar op cannot be coupled with non scalar op.')
    end
    
    
    

end

