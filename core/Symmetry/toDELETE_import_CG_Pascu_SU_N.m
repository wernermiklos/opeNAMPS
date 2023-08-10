function CGtensor= import_CG_Pascu_SU_N(filename)
%IMPORT_CG_Pascu_SU_N Summary of this function goes here
%   Detailed explanation goes here
    fileID = fopen(filename,'r');
    sizeRAW = [9, inf];
    RAW_data = fscanf(fileID,'%d %d %d %d %d %d %d %d %f',sizeRAW);
    sizeRAW = size(RAW_data);
    irrep_dims = containers.Map();
    outer_dims = containers.Map();
    blocks = containers.Map();
    
    for lineID = 1:sizeRAW(2)
        r1 = RAW_data(2,lineID) + 1;   % in MatLab we index from 1, not from 0
        m1 = RAW_data(3,lineID) + 1;
        r2 = RAW_data(4,lineID) + 1;
        m2 = RAW_data(5,lineID) + 1;
        R = RAW_data(6,lineID) + 1;
        M = RAW_data(7,lineID) + 1;
        alpha = RAW_data(8,lineID) + 1;
        key = char([r1,r2,R]);
        if ~isKey(irrep_dims,char(r1))
            irrep_dims(char(r1)) = m1;
        else
            irrep_dims(char(r1)) = max(m1,irrep_dims(char(r1)));
        end
        if ~isKey(irrep_dims,char(r2))
            irrep_dims(char(r2)) = m2;
        else
            irrep_dims(char(r2)) = max(m2,irrep_dims(char(r2)));
        end
        if ~isKey(irrep_dims,char(R))
            irrep_dims(char(R)) = M;
        else
            irrep_dims(char(R)) = max(M,irrep_dims(char(R)));
        end
        if ~isKey(outer_dims,key)
            outer_dims(key) = alpha;
        else
            outer_dims(key) = max(alpha,outer_dims(key));
        end
    end
    CGtensor = NAtensor({'m1', 'm2', 'M', 'alpha'}, {'i', 'i', 'o', 'o'} , ...
        {[1], [2], [3], [1, 2, 3]}, 1);
    
    
    for lineID = 1:sizeRAW(2)
        r1 = RAW_data(2,lineID) + 1;   % in MatLab we index from 1, not from 0
        m1 = RAW_data(3,lineID) + 1;
        r2 = RAW_data(4,lineID) + 1;
        m2 = RAW_data(5,lineID) + 1;
        R = RAW_data(6,lineID) + 1;
        M = RAW_data(7,lineID) + 1;
        alpha = RAW_data(8,lineID) + 1;
        value = RAW_data(9,lineID);
        key = char([r1,r2,R]);
        if ~isKey(blocks,key)
            blocks(key) = zeros(irrep_dims(key(1)),irrep_dims(key(2)),irrep_dims(key(3)),outer_dims(key));
        end
        tmp = blocks(key);    % I NEED TO DO THIS BECAUSE OF THE STUPID MAP CONTAINER!!!
        tmp(m1,m2,M,alpha) = value;
        blocks(key) = tmp;
    end
    
    for key_cell = keys(blocks)
        key = key_cell{1};
        CGtensor = NTset_block(CGtensor,{{'m1',1},{'m2',1},{'M',1}},{double(key(1)),double(key(2)),double(key(3))},{'m1','m2','M','alpha'},blocks(key));
    end
end

