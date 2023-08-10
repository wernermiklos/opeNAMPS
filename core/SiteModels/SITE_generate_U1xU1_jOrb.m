function out = SITE_generate_U1xU1_jOrb(rep,U1sym,info)
%GENERATE_JORB Summary of this function goes here
%   rep = 2J + 1

    if nargin == 2
        info = false;
    end
    eps = 1.e-12;
    jtot = (rep-1)/2;
    
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
    
    
    % We collect states into (rep_Q,rep_jz) sectors. sectors{rep_Q)(rep_jz) = [state1, state2, state3, ...]
    % rep_Q:  as usual the position in the series Q = 0, -1 ,+1, -2, +2,
    % ...
    % rep_jZ:    2*jz = 0, -1, +1, -2, +2, ...
    sectors = cell(1,rep+1);   % cell index represents charge, i.e. sectors{1} is for Q=0, ...
    num_of_states = cell(1,rep+1);
    for i = 1:(rep+1)
        num_of_states{i} = containers.Map('KeyType','double','ValueType','double');
        sectors{i} = containers.Map('KeyType','double','ValueType','any');
    end
    
    jzlist = jtot:-1:(-jtot);
    for key = statemap.keys
        jz = sum(jzlist(key{1}=='1'));
        Q = sum(key{1}=='1');
        rep_jz = Q_to_rep(2*jz);
        
        if ~isKey(sectors{Q+1},rep_jz)
            sectors{Q+1}(rep_jz) = [statemap(key{1})];
            num_of_states{Q+1}(rep_jz) = 1;
        else
            sectors{Q+1}(rep_jz) = [sectors{Q+1}(rep_jz),statemap(key{1})];
            num_of_states{Q+1}(rep_jz) = num_of_states{Q+1}(rep_jz) + 1;
        end
    end
    multiplet_list = {};
    
    for Q= 0:rep
        rep_Q = Q_to_rep(Q);
        for rep_jz_cell = num_of_states{Q+1}.keys
            rep_jz = rep_jz_cell{1};
            multiplet_list{end+1} = {[rep_Q, rep_jz], num_of_states{Q+1}(rep_jz)};
        end
    end

    out = SITE_create({U1sym,U1sym},multiplet_list,'FERMIONS',true);

    ph = NAtensor({'tau','mu','tau~','mu~'},{'o','o','i','i'},{[1],[1],[2],[2]},2);
    for secID = 1:length(multiplet_list)
        sector = multiplet_list{secID};
        multdim = 1;
        Q = rep_to_Q(sector{1}(1));
        ph = NTset_block(ph,{{'tau',1},{'tau~',1}},{sector{1},sector{1}},{'tau','mu','tau~','mu~'},reshape(((-1)^Q)*eye(sector{2}*multdim),[sector{2},multdim,sector{2},multdim]));
    end
    out = SITE_define_tensor_operator(out,'ph',{U1sym,U1sym},[1,1],{ph},1);
    if info
        disp('ph defined');
    end
    
    for i = 1:rep
        opjz = jzlist(i);
        rep_opjz = Q_to_rep(2*opjz);
        opQ = +1;
        rep_opQ = Q_to_rep(opQ);
        % cdag
        cdag_tmp = NAtensor({'tau','mu','tau~','mu~'},{'o','o','i','i'},{[1],[1],[2],[2]},2);
        for secID = 1:length(multiplet_list)
            sec = multiplet_list{secID};
            Q = rep_to_Q(sec{1}(1));
            jz = 0.5 * rep_to_Q(sec{1}(2));
            
            rep_jz = sec{1}(2);
            Qtild = Q + 1;
            jztild = jz + opjz;
            
            rep_Qtild = Q_to_rep(Qtild);
            rep_jztild = Q_to_rep(2*jztild);
     
            if SITE_sector_exist(out,[rep_Qtild,rep_jztild])
                sec_tild = {[rep_Qtild,rep_jztild],num_of_states{Qtild + 1}(rep_jztild)};
                cdag_block = full(cdag{i}(sectors{Qtild+1}(rep_jztild),sectors{Q+1}(rep_jz)));
                cdag_tmp = NTset_block(cdag_tmp,{{'tau',1},{'tau~',1}},...
                               {sec{1},sec_tild{1}},{'tau','mu','tau~','mu~'},...
                               reshape(cdag_block,[sec{2},1,sec_tild{2},1]));
            end
        end
        out = SITE_define_tensor_operator(out,['cdag_', num2str(rep_opjz)],{U1sym, U1sym},[rep_opQ,rep_opjz],{cdag_tmp},-1);    
    end
    
    for i = 1:rep
        opjz = - jzlist(i);
        rep_opjz = Q_to_rep(2*opjz);
        opQ = -1;
        rep_opQ = Q_to_rep(opQ);
        % cdag
        c_tmp = NAtensor({'tau','mu','tau~','mu~'},{'o','o','i','i'},{[1],[1],[2],[2]},2);
        for secID = 1:length(multiplet_list)
            sec = multiplet_list{secID};
            Q = rep_to_Q(sec{1}(1));
            jz = 0.5 * rep_to_Q(sec{1}(2));
            
            rep_jz = sec{1}(2);
            Qtild = Q + opQ;
            jztild = jz + opjz;
            
            rep_Qtild = Q_to_rep(Qtild);
            rep_jztild = Q_to_rep(2*jztild);
     
            if SITE_sector_exist(out,[rep_Qtild,rep_jztild])
                sec_tild = {[rep_Qtild,rep_jztild],num_of_states{Qtild + 1}(rep_jztild)};
                c_block = full(c{i}(sectors{Qtild+1}(rep_jztild),sectors{Q+1}(rep_jz)));
                c_tmp = NTset_block(c_tmp,{{'tau',1},{'tau~',1}},...
                               {sec{1},sec_tild{1}},{'tau','mu','tau~','mu~'},...
                               reshape(c_block,[sec{2},1,sec_tild{2},1]));
            end
        end
        replabel = Q_to_rep(-2*opjz);   % We label it to match the label of cdag
        out = SITE_define_tensor_operator(out,['c_', num2str(replabel)],{U1sym, U1sym},[rep_opQ,rep_opjz],{c_tmp},-1);
    end
    
end
    
function rep = Q_to_rep(Q)
    if abs(Q) < 0.00001
            rep = 1;
    else
            rep = round(2*abs(Q) + (Q + abs(Q))/(2*abs(Q)));
    end
end

function Q = rep_to_Q(rep)
    Q = (-1)^(rep-1) * floor(rep/2);
end
