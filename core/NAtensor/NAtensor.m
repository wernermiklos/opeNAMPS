function out = NAtensor(LegNames, LegTypes, Dependencies, NoOfSymmetries)
% Returns an NAtensor struct.
% ------------------
% LegNames:       {'name1', 'name2', ...}
% LegTypes:       {'type1', 'type2', ...}   type can be 'i' or 'o'
%                                           (incoming, outgoing)
% Dependencies:   {[dep1_1, dep1_2, ..], [dep2_1, dep2_2, ], ...}
% NoOfSymmetries: integer, the number of symmetries
           if ~iscell(LegNames) || ~iscell(LegTypes) || ~iscell(Dependencies)
              error(['Error in NAtensor: LegNames, LegTypes and Dependencies', ...
                  'must be cell arrays'])
           end
           
           if NoOfSymmetries < 1
              error('Error in NAtensor: NoOfSymmetries must be a positive integer')
           end
           
           if length(unique(LegNames)) ~= length(LegNames)
              error('Error in NAtensor: redundancy in LegNames')
           end
           
           if length(LegNames) ~= length(LegTypes)
              error('Error in NAtensor: length of LegTypes differs from length of LegNames')
           end
           
           if length(LegNames) ~= length(Dependencies)
              error('Error in NAtensor: length of Dependencies differs from length of LegNames')
           end
           
           for legtype = LegTypes
               if legtype{1} ~= 'i' && legtype{1} ~= 'o'
                   error('Error in NAtensor: legtypes can only be "i" or "o".')
               end
           end
           
           test_dep_set = Dependencies{1};
           for i = 2:length(LegNames)
               test_dep_set = union(test_dep_set,Dependencies{i});
           end
           if size(test_dep_set,1) > 1     %we make sure this is a row vector (union() returns a column vector if one of the sets is empty...)
               test_dep_set = test_dep_set';
           end
           
           if ~isempty(test_dep_set)
               if ~isequal(test_dep_set,1:test_dep_set(length(test_dep_set)))
                   error(['Error in NAtensor: irrep indices are not the integers 1,2,3,4,5,...! ' , ...
                       mat2str(test_dep_set)])
               end
           end
           
           out.type = 'NAtensor';
           out.leg_names = LegNames;
           out.leg_types = LegTypes;
           out.data = containers.Map('UniformValues',false);
           out.no_of_symmetries = NoOfSymmetries;
           out.irrep_number = length(test_dep_set);
           out.dependencies = Dependencies;

           % out.shape = containers.Map('UniformValues',false);  with the
           %                         new (c-mex) tensordot we don't need it anyore. 
end

