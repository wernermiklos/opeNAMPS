function Coupling_RED = COUP_extract_RED(Coupling_FULL, CG, PANDA, center_pos,no_of_sites,varargin)
%COUP_extract_RED Summary of this function goes here
%   Coupling_FULL:      The full coupling tensor (with "M" and "eta" legs)
%   CG:                 The (total) Clebsch-tensor of the symmetry (can be an LNAtensor for no_of_symmetries>1)
%   PANDA:              The (total) PANDA-tensor (projection of CG) of the symmetry.
%   center_pos:         the center position (the "the last leg of the left part")


    parameters = struct();
    parameters.OPS_ARE_SCALAR = false;
    parameters = parameter_updater(parameters,varargin);

    eta_groups = get_eta_legs(Coupling_FULL,no_of_sites);
    new_eta_names = cell(1,length(eta_groups));

    j = 0;
    ops_are_scalar = parameters.OPS_ARE_SCALAR;

    if center_pos < 1 || center_pos > no_of_sites
        error('center_pos must be between 1 and no_of_sites')
    end


    for i = 1:(center_pos)
        if ~isempty(eta_groups{i})
            j = j + 1;
            new_eta_names{i} = cell(1,length(eta_groups{i}));
            for k = 1:length(eta_groups{i})
                new_eta_names{i}{k} = ['L_', eta_groups{i}{k}(1:(end-1)), num2str(j)];
            end
            Coupling_FULL = NTrename_legs(Coupling_FULL,eta_groups{i},new_eta_names{i});
        end
    end

    j = 0;
    for i = no_of_sites:-1:(center_pos+1)
        if ~isempty(eta_groups{i})
            j = j + 1;
            new_eta_names{i} = cell(1,length(eta_groups{i}));
            for k = 1:length(eta_groups{i})
                new_eta_names{i}{k} = ['R_', eta_groups{i}{k}(1:(end-1)), num2str(j)];
            end
            Coupling_FULL = NTrename_legs(Coupling_FULL,eta_groups{i},new_eta_names{i});
        end
    end


    Coupling_RED = Coupling_FULL;
    if ~ops_are_scalar
        if center_pos == no_of_sites || center_pos == no_of_sites-1
            center_pos = no_of_sites-1;   % From here the two cases are equivalent
        end
        Coupling_RED = NTrename_legs(Coupling_RED,{'M_1'},{'L_M'});
        Coupling_RED = NTrename_legs(Coupling_RED,{['M_',num2str(no_of_sites)]},{'R_M'});
        for i = 2:center_pos
            Coupling_RED = NTdot(Coupling_RED,...
                                 NTconj(CG),...
                                 {'L_M',['M_',num2str(i)]},...
                                 {'m1','m2'},...
                                 {{},{{'alpha',['L_xi_',num2str(i-1)]},{'M','L_M'}}});
        end
        for i = (no_of_sites-1):-1:(center_pos+1)
            Coupling_RED = NTdot(Coupling_RED,...
                                 NTconj(CG),...
                                 {'R_M',['M_',num2str(i)]},...
                                 {'m1','m2'},...
                                 {{},{{'alpha',['R_xi_',num2str(no_of_sites-i)]},{'M','R_M'}}});
        end
        Coupling_RED = NTdot(Coupling_RED,...
                             NTconj(PANDA),...
                             {'L_M','R_M'},...
                             {'m1','m2'});
    end
end



% ---------------- Auxiliary functions -----------------

function outleggroups = get_eta_legs(Coupling_FULL,no_of_sites)
    outleggroups = cell(1,no_of_sites);
    for name = Coupling_FULL.leg_names
        if length(name{1}) > 2
            if strcmp(name{1}(1:3),'eta')
                n = str2num(name{1}(end));
                if isempty(outleggroups{n})
                    outleggroups{n} = {name{1}};
                else
                    outleggroups{n}{end+1} = name{1};
                end
            end
        end
    end
end