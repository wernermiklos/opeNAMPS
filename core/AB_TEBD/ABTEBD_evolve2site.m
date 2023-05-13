function [ABTEBD_env, info] = ABTEBD_evolve2site(ABTEBD_env, bond_pos, evolver,Mmax, varargin)
%ABTEBD_APPLYEVOLVER Summary of this function goes here
%   Detailed explanation goes here
parameters = struct();
parameters.NORMALIZE = true;
parameters.TRUNCATE_AT_MULTIPLETS=false;
parameters = parameter_updater(parameters,varargin);


A1 = NTrename_legs(ABTEBD_env.StateMPS.left_matrices{bond_pos},...
                   {'t_in','tau','t_out'},{'t_left','tau_1','t_right'});
A2 = NTrename_legs(ABTEBD_env.StateMPS.left_matrices{mod(bond_pos,ABTEBD_env.chain_length)+1}, ...
                   {'t_in','tau','t_out'},{'t_left','tau_2','t_right'});
l2 = ABTEBD_env.StateMPS.schmidt_list{mod(bond_pos,ABTEBD_env.chain_length)+1};

W = NTdot(A1,A2,{'t_right'},{'t_left'});
Wnew = NTdot(W,evolver,{'tau_1','tau_2'},{'tau_1','tau_2'},{{},{{'tau_1~','tau_1'},{'tau_2~','tau_2'}}});
Theta_new = NTdot(Wnew,l2,{'t_right'},{'t_left'});
if parameters.NORMALIZE
    normsquare = NTdot(Theta_new, NTconj(Theta_new), {'t_left', 'tau_1', 'tau_2', 't_right'},  {'t_left', 'tau_1', 'tau_2', 't_right'});
    Theta_new = NTmult(Theta_new,1.0/sqrt(normsquare));
    Wnew = NTmult(Wnew,1.0/sqrt(normsquare));
end

svd_output = ABTEBD_svd(Theta_new,ABTEBD_env.Symmetries,Mmax);



info.normsq = svd_output.normsq;
info.normsq_untr = svd_output.normsq_untr;
info.trunc_err = 1.0 - info.normsq/info.normsq_untr;
info.bond_sectors = NTget_leg_sectors(svd_output.schmidt,'t_left');
info.bond_dim = 0;
for i = 1:length(info.bond_sectors)
    info.bond_dim = info.bond_dim + info.bond_sectors{i}{2};
end



if parameters.TRUNCATE_AT_MULTIPLETS && info.bond_dim==Mmax
    [irr,val]=NTexport_all_data_blocks(svd_output.schmidt,{{'t_left',1},{'t_right',1}},{'t_left','t_right'});
    schmidt_val=[];
    for ind=1:length(irr)
        schmidt_val(end+1:end+length(val{ind}))=diag(val{ind});
    end
    m=min(schmidt_val);

    schm=NAtensor({'t_left','t_right'},{'i','i'},{[1],[2]},svd_output.schmidt.no_of_symmetries);
    id_cut_left=NAtensor({'t_in','t_out'},{'i','o'},{[1],[2]},svd_output.schmidt.no_of_symmetries);
    for ind=1:length(irr)
        indmx=double(((val{ind}-m)>10^(-8) | (val{ind}-m)/m > 0.01));
       schm = NTset_block(schm,{{'t_left',1},{'t_right',1}},{irr{ind}{1},irr{ind}{2}},{'t_left','t_right'},val{ind}.*indmx);
       id_cut_left = NTset_block(id_cut_left,{{'t_in',1},{'t_out',1}},{irr{ind}{1},irr{ind}{1}},{'t_in','t_out'},indmx);
    end
    svd_output.A = NTdot(svd_output.A,id_cut_left,{'t_out'},{'t_in'});      
    svd_output.schmidt = schm;
    svd_output.normsq = NTdot(schm,NTconj(schm),{'t_left','t_right'},{'t_left','t_right'});
end



norm_factor = sqrt(svd_output.normsq_untr / svd_output.normsq);
ABTEBD_env.StateMPS.left_matrices{mod(bond_pos,ABTEBD_env.chain_length)+1} = ... 
                NTmult(NTdot(NTconj(svd_output.A),Wnew, ...
                     {'t_in','tau'},...
                     {'t_left','tau_1'},...
                     {{{'t_out','t_in'}},...
                      {{'t_right','t_out'},{'tau_2','tau'}}}), ...
               norm_factor);
ABTEBD_env.StateMPS.left_matrices{bond_pos} = svd_output.A;
ABTEBD_env.StateMPS.schmidt_list{bond_pos} = NTmult(svd_output.schmidt, norm_factor);
if bond_pos == ABTEBD_env.chain_length
            ABTEBD_env.StateMPS = ABMPS_set_schmidt(ABTEBD_env.StateMPS,ABTEBD_env.StateMPS.schmidt_list{bond_pos});
end


end


