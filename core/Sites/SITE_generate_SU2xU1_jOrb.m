function out = SITE_generate_SU2xU1_jOrb(rep,SU2sym,U1sym,info)
%GENERATE_JORB Summary of this function goes here
%   rep = 2J + 1

    if nargin == 3
        info = false;
    end
    eps = 1.e-12;
    
    statemap = containers.Map('KeyType','char','ValueType','double');  %Will translate binary string like '0100011101' to state index. 
    
    
    for i = 0:(2^rep-1)
        key = dec2bin(i);    % binary key (occupation number rep) '010100110101'
        tmp = length(key);
        for j = 1:(rep-tmp)
            key = ['0',key];   % We extend the binary string with zeros from the left, if necessary
        end
        statemap(key) = i+1;
    end
    
    
    % We build up the cdag operator
    cdag = cell(1,rep); % cells for the 2j+1 f_dagger operators;
    c = cell(1,rep); % cells for the 2j+1 f operators;
    for i = 1:rep
        cdag{i} = sparse(2^rep,2^rep);
        for key = statemap.keys
            if key{1}(i) ~= '1'
                val = 1;
                for j = 1:(i-1)
                    if key{1}(j) == '1'
                        val = -val;
                    end
                end
                outkey = key{1};
                outkey(i) = '1';
                cdag{i}(statemap(outkey),statemap(key{1})) = val;
            end
        end
        c{i} = cdag{i}';
    end

    
    % We build up the spin operators
    Sz = sparse(2^rep,2^rep);
    Sp = sparse(2^rep,2^rep);
    Sm = sparse(2^rep,2^rep);
    jtot = (rep-1)/2;
    for i = 1:rep
        jz= jtot-(i-1);
        Sz = Sz + jz*cdag{i}*(c{i});
        if i ~= rep
            Sm = Sm + sqrt(jtot*(jtot+1)-jz*(jz-1))*cdag{i+1}*(c{i});
        end
        if i ~= 1
            Sp = Sp + sqrt(jtot*(jtot+1)-jz*(jz+1))*cdag{i-1}*(c{i});
        end
    end
    Ssq = Sz*Sz + 0.5*(Sp*Sm + Sm*Sp);
    
    
    % We collect states into (Q,Jz) sectors. sectors{Q+1)(Jz) = [state1, state2, state3, ...]
    sectors = cell(1,rep+1);   % cell index represents charge, i.e. sectors{1} is for Q=0, ...
    for i = 1:(rep+1)
        sectors{i} = containers.Map('KeyType','double','ValueType','any');
    end
    num_of_states = zeros(1,rep+1);
    jzlist = jtot:-1:(-jtot);
    for key = statemap.keys
        jz = sum(jzlist(key{1}=='1'));
        Q = sum(key{1}=='1');
        state = zeros(2^rep,1);
        num_of_states(Q+1) = num_of_states(Q+1)+1;
        state(statemap(key{1}))=1;
        if ~isKey(sectors{Q+1},jz)
            sectors{Q+1}(jz) = state;
        else
            sectors{Q+1}(jz) = [sectors{Q+1}(jz),state];
        end
    end
    
    % Finally we build up the multiplets
    multiplets = cell(1,rep+1);   % cell index is Q+1, i.e. sectors{1} stands for Q=0, ...
    for i = 1:(rep+1)
        multiplets{i} = containers.Map('KeyType','double','ValueType','any');
    end
    
    U1replist = 1:2:(2*rep+1);
    SU2replist = [];
    
    for Q = 0:rep
        remaining = num_of_states(Q+1);   % The number of remaining states in the Q sector.
        while remaining > 0
            jzmax = max(cell2mat(sectors{Q+1}.keys));    % The maximal jz value found
            jzmax_states = sectors{Q+1}(jzmax);          % States with maximal jz value
            startstate = jzmax_states(:,1);              % We start the iteration from this state
            sectors{Q+1}(jzmax) = jzmax_states(:,2:end);     % we remove the startstate from the pool.
            remaining = remaining - 1;
            
            shape = size(sectors{Q+1}(jzmax));               
            if shape(2) == 0    % If the (Q, jzmax) sector is empty, we delete it.
                sectors{Q+1}.remove(jzmax);
            end
            
            multiplet = zeros(2^rep,2*jzmax+1);     % We allocate space for the multiplet
            multiplet(:,1) = startstate;            % We put startstate into the first column.
            for i = 2:(2*jzmax+1)
                newstate = Sm*multiplet(:,i-1);     % We apply S_minus operator 
                newstate = newstate/norm(newstate);  
                samejz_states = sectors{Q+1}(jzmax-i+1);  % We take the states with the same jz as newstate;
                overlaps = newstate'*samejz_states;
                
                
                samejz_states = samejz_states - newstate*overlaps;      % We remove the parallel parts from the vectors                

%               ALTERNATIVE ORTHOGONALIZATION:                
%                 M = samejz_states'*samejz_states;
%                 [v,w] = eig(M);
%                 newjz_states = samejz_states*v;
%                 for state_id = 1:size(newjz_states,2)
%                     if norm(newjz_states(:,state_id)) > eps
%                        newjz_states(:,state_id) = newjz_states(:,state_id)/norm(newjz_states(:,state_id));
%                     end
%                 end
% 
%                 S = diag(w);
                [newjz_states, S] = svd(samejz_states,'econ');          % We orthogonalize using svd.
                shape = size(newjz_states);
                 if S(shape(2),shape(2)) > eps                           % If the discarded Schmidt weight is larger than eps, we raise a Warning
                     warning('Warning, SVD')
                     S(shape(2),shape(2))
                 end
                % newjz_states = gram_schmidt_orth([newstate, samejz_states]);   % MatLab built in "orth" seems to be better...   
                
                
                newjz_states = newjz_states(:,1:(end-1));                     % We discard "newstate")
                sectors{Q+1}(jzmax-i+1) = newjz_states;
                shape = size(sectors{Q+1}(jzmax-i+1));
                if shape(2) == 0    % If the (Q, jz) sector is empty, we delete it.
                    sectors{Q+1}.remove(jzmax-i+1);
                end
                remaining = remaining - 1;
                multiplet(:,i) = newstate;
            end
            
            % Finally we store the new multiplet. Multiplets are indexed by
            % 2*jzmax + 1 (that is the "representation" index.)
            if ~isKey(multiplets{Q+1},2*jzmax+1)
                multiplets{Q+1}(2*jzmax+1) = {multiplet};
            else
                multiplets{Q+1}(2*jzmax+1) = [multiplets{Q+1}(2*jzmax+1),{multiplet}]; 
            end
            SU2replist = [SU2replist,2*jzmax+1];
        end    
    end
    
    SU2replist = unique(SU2replist);
   
    %CG_SU2 = generate_SU2_CG(SU2replist,1:(floor(rep/2)*max(SU2replist)),SU2replist);
    %CG_U1 = generate_U1_CG(U1replist,1:max(U1replist),U1replist);
    SU2replist_test1 = NTget_irrep_set(SU2sym.CGtensor,{'m1',1});
    SU2replist_test2 = NTget_irrep_set(SU2sym.CGtensor,{'m2',1});
    SU2replist_test3 = NTget_irrep_set(SU2sym.CGtensor,{'M',1});
    U1replist_test1 = NTget_irrep_set(U1sym.CGtensor,{'m1',1});
    U1replist_test2 = NTget_irrep_set(U1sym.CGtensor,{'m2',1});
    U1replist_test3 = NTget_irrep_set(U1sym.CGtensor,{'M',1});
    for SU2rep = SU2replist
        found1 = any(cellfun(@(x) eq(SU2rep,x),SU2replist_test1));
        found3 = any(cellfun(@(x) eq(SU2rep,x),SU2replist_test3));
        if ~found1 || ~found3
            error(['SU2sym does not contain rep ', mat2str(SU2rep), ' on leg m1 or M'])
        end
    end
    for SU2rep = 1:(2*rep-1)                 % Will be enough for 2-particle interaction  (max operator rep on size is 2*j)
        found2 = any(cellfun(@(x) eq(SU2rep,x),SU2replist_test2));
        if ~found2
            error(['SU2sym does not contain rep ', mat2str(SU2rep), ' on leg m2'])
        end
    end
    for U1rep = 1:2:(2*rep+1)
        found1 = any(cellfun(@(x) eq(U1rep,x),U1replist_test1));
        found3 = any(cellfun(@(x) eq(U1rep,x),U1replist_test3));
        if ~found1 || ~found3
            error(['U1sym does not contain rep ', mat2str(U1rep), ' on leg m1 or M'])
        end
    end
    for U1rep = 1:max(U1replist)
        found2 = any(cellfun(@(x) eq(U1rep,x),U1replist_test2));
        if ~found2
            error(['U1sym does not contain rep ', mat2str(U1rep), ' on leg m2'])
        end
    end

   
    
    multiplet_list = {};
    for Q = 0:rep
        for SU2rep_cell = multiplets{Q+1}.keys
            SU2rep = SU2rep_cell{1};
            U1rep = 2*Q+1;
            multiplet_list{end+1} = {[SU2rep,U1rep],length(multiplets{Q+1}(SU2rep))};
        end
    end
    
    out = SITE_create({SU2sym,U1sym},multiplet_list,'FERMIONS',true);
    ph = NAtensor({'tau','mu','tau~','mu~'},{'o','o','i','i'},{[1],[1],[2],[2]},2);
    for secID = 1:length(multiplet_list)
        sector = multiplet_list{secID};
        multdim = SU2sym.irrep_dimensions(sector{1}(1));
        Q = (sector{1}(2)-1)/2;
        ph = NTset_block(ph,{{'tau',1},{'tau~',1}},{sector{1},sector{1}},{'tau','mu','tau~','mu~'},reshape(((-1)^Q)*eye(sector{2}*multdim),[sector{2},multdim,sector{2},multdim]));
    end
    fdag_list = cell(1,rep);
    f_list = cell(1,rep);
    f_signs = NTget_block(SU2sym.CGtensor,{{'m1',1},{'m2',1},{'M',1}},{rep,rep,1},{'m1','m2','M','alpha'});
    f_signs = sign(f_signs);
    for M_op = 1:rep
        fdag_list{M_op} = NAtensor({'tau','mu','tau~','mu~'},{'o','o','i','i'},{[1],[1],[2],[2]},2);
        f_list{M_op} = NAtensor({'tau','mu','tau~','mu~'},{'o','o','i','i'},{[1],[1],[2],[2]},2);
        f_sign = f_signs(M_op,rep-M_op+1)*f_signs(1,rep);    % Sign convention f{1} = + c{rep}
        for Q = 1:rep
            U1rep1 = 2*Q-1;
            U1rep2 = 2*Q+1;
            for SU2rep1_cell = multiplets{Q}.keys
                SU2rep1 = SU2rep1_cell{1};
                for SU2rep2_cell = multiplets{Q+1}.keys
                    SU2rep2 = SU2rep2_cell{1};
                    multiplet_list1 = multiplets{Q}(SU2rep1);
                    multiplet_list2 = multiplets{Q+1}(SU2rep2);
                    irrepdim1 = out.multiplet_dims(char([SU2rep1, U1rep1]));
                    irrepdim2 = out.multiplet_dims(char([SU2rep2, U1rep2]));
                    
                    newblock_fdag = zeros(length(multiplet_list1),irrepdim1,length(multiplet_list2),irrepdim2);
                    newblock_f = zeros(length(multiplet_list2),irrepdim2,length(multiplet_list1),irrepdim1);
                    for tau = 1:length(multiplet_list1)
                        for tautild = 1:length(multiplet_list2)
                            newblock_fdag(tau,:,tautild,:) = ((multiplet_list2{tautild}')*cdag{M_op}*multiplet_list1{tau}).';
                            newblock_f(tautild,:,tau,:) = f_sign*((multiplet_list1{tau}')*c{rep-M_op+1}*multiplet_list2{tautild}).';
                        end
                    end
                    fdag_list{M_op} = NTset_block(fdag_list{M_op}, {{'tau',1},{'tau~',1}},{[SU2rep1,U1rep1],[SU2rep2,U1rep2]},{'tau','mu','tau~','mu~'},newblock_fdag);
                    f_list{M_op} = NTset_block(f_list{M_op}, {{'tau',1},{'tau~',1}},{[SU2rep2,U1rep2],[SU2rep1,U1rep1]},{'tau','mu','tau~','mu~'},newblock_f);
                end
            end
        end
    end
    
    opirrep_fdag = [rep, 3];    %rep = 2j + 1,   U1rep = 3 means charge = +1
    opirrep_f = [rep, 2];    %rep = 2j + 1,   U1rep = 2 means charge = -1
    out = SITE_define_tensor_operator(out, 'fdag',{SU2sym,U1sym},opirrep_fdag,fdag_list,'FERMIONICITY',-1);
    if info
        disp('fdag defined');
    end
    out = SITE_define_tensor_operator(out,'f',{SU2sym,U1sym},opirrep_f,f_list,'FERMIONICITY',-1);
    if info
        disp('f defined');
    end
    out = SITE_define_tensor_operator(out,'ph',{SU2sym,U1sym},[1,1],{ph},'FERMIONICITY',1);
    if info
        disp('ph defined');
    end
    out = SITE_dot_define(out, 'fdag','fdag','fdag_fdag',{SU2sym,U1sym});
    if info
        disp('fdag_fdag defined');
    end
    out = SITE_dot_define(out,'fdag','f','fdag_f',{SU2sym,U1sym});
    if info
        disp('fdag_f defined');
    end
    out = SITE_dot_define(out,'f','f','f_f',{SU2sym,U1sym});
    if info
        disp('f_f defined');
    end
    out = SITE_dot_define(out,'fdag_fdag', 'f', 'fdag_fdag_f', {SU2sym,U1sym});
    if info
        disp('fdag_fdag_f defined')
    end
    out = SITE_dot_define(out,'fdag_f', 'f', 'fdag_f_f', {SU2sym,U1sym});
    if info
        disp('fdag_f_f defined')
    end
    out = SITE_dot_define(out,'fdag_fdag_f', 'f', 'fdag_fdag_f_f', {SU2sym,U1sym},{[1,1]});   % From this product we keep only the trivial [1,1] scalar operator.
    if info
        disp('fdag_fdag_f_f defined')
    end
end
    
    
    
