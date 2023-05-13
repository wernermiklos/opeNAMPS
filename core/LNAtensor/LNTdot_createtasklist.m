function [keylist, tasklist] = LNTdot_createtasklist(objkey_lists, ...
                                                     otherkeys, ...
                                                     matchkey_layout_list1, ...
                                                     matchkey_layout_list2, ...
                                                     outkey_layout_list1, ...
                                                     outkey_layout2, ...
                                                     outkey_length, ...
                                                     no_of_symmetries)
%LNTDOT_CREATETASKLIST Summary of this function goes here
%   Detailed explanation goes here
    if ~verLessThan('matlab','9.13')
        matchkey_table1 = cellfun(@(x) dictionary(strings(0),[]), cell(1, no_of_symmetries), 'UniformOutput', false);
        bucket_list = cell(1,no_of_symmetries);
        for symID = 1:no_of_symmetries
            objkey_tmplist = objkey_lists{symID};
            matchkey_layout_tmp = matchkey_layout_list1{symID};
            bucket_list{symID} = cell(1,length(objkey_tmplist));
            last_bucket = 0;
            for i = 1:length(objkey_tmplist)
                objmatchkey = objkey_tmplist{i}(matchkey_layout_tmp);
                if ~isKey(matchkey_table1{symID},objmatchkey)
                    last_bucket = last_bucket + 1;
                    matchkey_table1{symID}(objmatchkey) = last_bucket;
                    bucket_list{symID}{last_bucket} = i;
                else
                    bucketID = matchkey_table1{symID}(objmatchkey);
                    bucket_list{symID}{bucketID} = [bucket_list{symID}{bucketID},i];
                end
            end
        end
        
        outkey_map = dictionary(strings(0),[]);
        tasklist = zeros(0,no_of_symmetries+2);
        outcount = 0;
        for bID2 = 1:length(otherkeys)
            nonzero = true;
            otherkey = otherkeys{bID2};
            outkey = char(zeros(1,outkey_length));
            outkey(outkey_layout2(outkey_layout2~=0)) = otherkey(outkey_layout2 ~= 0);   % We set the known part of outkey (that comes from otherkey)
            obj_bIDs = cell(1,no_of_symmetries);
            for symID = 1:no_of_symmetries
                matchkey2 = otherkey(matchkey_layout_list2{symID});
                if isKey(matchkey_table1{symID}, matchkey2)
                    bucketID = matchkey_table1{symID}(matchkey2);
                    obj_bIDs{symID} = bucket_list{symID}{bucketID};
                else
                    nonzero = false;
                    break;
                end
            end
            if ~nonzero
                continue
            end
            bIDcombs = cartprod_fromcell(obj_bIDs)';  
            tasklist_tmp = [bIDcombs,bID2*ones(size(bIDcombs,1),1),zeros(size(bIDcombs,1),1)];
            sub_bIDs = zeros(1,no_of_symmetries);
            for taskID = 1:size(tasklist_tmp,1)
                for symID = 1:no_of_symmetries
                    if tasklist_tmp(taskID,symID) ~= sub_bIDs(symID);
                        sub_bIDs(symID) = tasklist_tmp(taskID,symID);
                        subkey = objkey_lists{symID}{sub_bIDs(symID)};
                        outkey(outkey_layout_list1{symID}(outkey_layout_list1{symID}~=0)) = subkey(outkey_layout_list1{symID}~=0);
                    end
                end
                if ~isKey(outkey_map,outkey)
                    outcount = outcount + 1;
                    outkey_map(outkey) = outcount;
                    tasklist_tmp(taskID,end) = outcount;
                else
                    tasklist_tmp(taskID,end) = outkey_map(outkey);
                end
            end
            tasklist = [tasklist; tasklist_tmp];
        end
        keylist = strings(1,outcount);
        outkeys = keys(outkey_map);
        outindices = values(outkey_map);
        keylist(outindices) = outkeys;
        keylist = cellfun(@(x) char(x), num2cell(keylist), 'UniformOutput',false);
        tasklist = sortrows(tasklist);
    else
        matchkey_table1 = cellfun(@(x) containers.Map('UniformValues',false), cell(1, no_of_symmetries), 'UniformOutput', false);
        for symID = 1:no_of_symmetries
            objkey_tmplist = objkey_lists{symID};
            matchkey_layout_tmp = matchkey_layout_list1{symID};
            for i = 1:length(objkey_tmplist)
                objmatchkey = objkey_tmplist{i}(matchkey_layout_tmp);
                if ~isKey(matchkey_table1{symID},objmatchkey)
                    matchkey_table1{symID}(objmatchkey) = i;
                else
                    matchkey_table1{symID}(objmatchkey) = [matchkey_table1{symID}(objmatchkey),i];
                end
            end
        end
        
        outkey_map = containers.Map('KeyType','char','ValueType','double');
        tasklist = zeros(0,no_of_symmetries+2);
        outcount = 0;
        for bID2 = 1:length(otherkeys)
            nonzero = true;
            otherkey = otherkeys{bID2};
            outkey = char(zeros(1,outkey_length));
            outkey(outkey_layout2(outkey_layout2~=0)) = otherkey(outkey_layout2 ~= 0);   % We set the known part of outkey (that comes from otherkey)
            obj_bIDs = cell(1,no_of_symmetries);
            for symID = 1:no_of_symmetries
                matchkey2 = otherkey(matchkey_layout_list2{symID});
                if isKey(matchkey_table1{symID}, matchkey2)
                    obj_bIDs{symID} = matchkey_table1{symID}(matchkey2);
                else
                    nonzero = false;
                    break;
                end
            end
            if ~nonzero
                continue
            end
            bIDcombs = cartprod_fromcell(obj_bIDs)';  
            tasklist_tmp = [bIDcombs,bID2*ones(size(bIDcombs,1),1),zeros(size(bIDcombs,1),1)];
            sub_bIDs = zeros(1,no_of_symmetries);
            for taskID = 1:size(tasklist_tmp,1)
                for symID = 1:no_of_symmetries
                    if tasklist_tmp(taskID,symID) ~= sub_bIDs(symID);
                        sub_bIDs(symID) = tasklist_tmp(taskID,symID);
                        subkey = objkey_lists{symID}{sub_bIDs(symID)};
                        outkey(outkey_layout_list1{symID}(outkey_layout_list1{symID}~=0)) = subkey(outkey_layout_list1{symID}~=0);
                    end
                end
                if ~isKey(outkey_map,outkey)
                    outcount = outcount + 1;
                    outkey_map(outkey) = outcount;
                    tasklist_tmp(taskID,end) = outcount;
                else
                    tasklist_tmp(taskID,end) = outkey_map(outkey);
                end
            end
            tasklist = [tasklist; tasklist_tmp];
        end
        outkeys = keys(outkey_map);
        outindices = cell2mat(values(outkey_map));
        keylist(outindices) = outkeys;
        tasklist = sortrows(tasklist);
    end
end

