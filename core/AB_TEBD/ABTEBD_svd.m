function outstruct = ABTEBD_svd(psi, symmetries, Mmax)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    outstruct = struct();
    outstruct.normsq_untr = NTdot(psi,NTconj(psi), ...
                                  {'t_left','tau_1','tau_2','t_right'}, ...
                                  {'t_left','tau_1','tau_2','t_right'});
    svdten = NAtensor({'t_left','tau_1','tau_2', 't_right'},{'i','i','i','i'},{[1,5],[2],[3],[4,5]},psi.no_of_symmetries);

    [irreps, blocks,~] = NTexport_all_data_blocks(psi,{{'t_left',1},{'tau_1',1},{'tau_2',1},{'t_right',1}},{'t_left','tau_1','tau_2', 't_right'});
    newkeys = cell(1,length(blocks));
    
    %we add formally the 5th irrep (the center irrep) as a dependency if
    %the t_left and t_right legs
    for bID = 1:length(blocks)
        Gamma_left = irreps{bID}{1};
        Gamma_1 = irreps{bID}{2};
        Gamma_2 = irreps{bID}{3};
        Gamma_right = irreps{bID}{4};
        fusion = SUPER_fusion_rule(symmetries,Gamma_left,Gamma_1);
        Gamma_center = fusion{1}{1};
        newkeys{bID} = char([Gamma_left,Gamma_1,Gamma_2,Gamma_right,Gamma_center]);
    end
    if ~isempty(blocks)
        svdten.data = containers.Map(newkeys,blocks,'UniformValues',false);
    end
    [U,newschmidt,V,normsq] = NTsvd(svdten,{'t_left', 'tau_1'},{'tau_2', 't_right'},{'t_bond_left', 't_bond_right'},{'i','i'},Mmax);
    
    %here we delete first the dummy (auxiliary) dependencies back
    U.dependencies{find(strcmp(U.leg_names,'t_left'),1)}(2) = [];
    V.dependencies{find(strcmp(V.leg_names,'t_right'),1)}(2) = [];
    
    %we switch the irrep of t_bond_right to the conjugate rep in newschmidt
    oldkeys = keys(newschmidt.data);
    blocks = values(newschmidt.data);
    newkeys = cell(1,length(blocks));
    for bID = 1:length(blocks)
         newkeys{bID} = char([oldkeys{bID},SUPER_conj_rep(symmetries,oldkeys{bID})]);
    end
    newschmidt.irrep_number = 2;
    newschmidt.dependencies{find(strcmp(newschmidt.leg_names,'t_bond_right'),1)} = [2];
    if ~isempty(blocks)
        newschmidt.data = containers.Map(newkeys,blocks,'UniformValues',false);
    end
    
    %we switch the irrep of t_bond_right (the 3rd irrep in V) to the conjugate rep in V
    oldkeys = keys(V.data);
    blocks = values(V.data);
    newkeys = cell(1,length(blocks));
    for bID = 1:length(blocks)
         newkeys{bID} = char([oldkeys{bID}(1:(2*V.no_of_symmetries)), ...
                              SUPER_conj_rep(symmetries,oldkeys{bID}((2*V.no_of_symmetries+1):end))]);
    end
    if ~isempty(blocks)
        V.data = containers.Map(newkeys,blocks,'UniformValues',false);
    end
    
    outstruct.schmidt = NTrename_legs(newschmidt, {'t_bond_left', 't_bond_right'},{'t_left','t_right'});
    outstruct.A = NTrename_legs(U,{'t_left','tau_1','t_bond_left'},{'t_in','tau','t_out'});
    outstruct.B = NTrename_legs(V,{'t_bond_right','tau_2','t_right'},{'t_out','tau','t_in'});
    outstruct.normsq = normsq;
          


end

