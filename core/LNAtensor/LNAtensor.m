function obj = LNAtensor(LegNames, LegTypes, Dependencies, NoOfSymmetries)
           % Constructor of LNAtensor objects. 
           % (Works similarly to standard NAtensors)
           % ------------------
           % obj = LNAtensor(LegNames, LegTypes, Dependencies, NoOfSymmetries)
           % ------------------
           % LegNames:       {'name1', 'name2', ...}
           % LegTypes:       {'type1', 'type2', ...}   type can be 'i' or 'o'
           % Dependencies:   {[dep1_1, dep1_2, ..], [dep2_1, dep2_2, ], ...}
           % NoOfSymmetries: integer, the number of symmetries
           if ~iscell(LegNames) || ~iscell(LegTypes) || ~iscell(Dependencies)
              error(['Error in NAtensor: LegNames, LegTypes and Dependencies', ...
                  'must be cell arrays'])
           end
           
           if NoOfSymmetries < 1
              error('Error in NAtensor: NoOfSymmetries must be a positive integer')
           elseif NoOfSymmetries == 1
               warning('Instead of LNAtensor an NAtensor is returned if NoOfSymmetries = 1')
               obj = NAtensor(LegNames,LegTypes,Dependencies,NoOfSymmetries);
               return;
           else
               obj.sub_tensors = cell(1,NoOfSymmetries);
               for symID = 1:NoOfSymmetries
                   obj.sub_tensors{symID} = NAtensor(LegNames, LegTypes, Dependencies, 1);   % no_of_symmetries = 1 for each subtensor
               end
               obj.type = 'LNAtensor';
               obj.leg_names = LegNames;
               obj.leg_types = LegTypes;
               obj.no_of_symmetries = NoOfSymmetries;
               obj.irrep_number = obj.sub_tensors{1}.irrep_number;
               obj.dependencies = Dependencies;
               obj.sub_just_ones = true(1,NoOfSymmetries);  % Initially this is true (empty sub tensors)
               obj.sub_just_nums = true(1,NoOfSymmetries);  % Initially this is true (empty sub tensors)
           end
end


