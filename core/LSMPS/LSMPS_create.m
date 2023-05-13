function obj = LSMPS_create(chain_length, no_of_symmetries)
   % Constructor of LeftSymmetric MPS.
   % LSMPS uses the convention where the Clebsch layer is always fully left
   % canonical, and bond indices are characterized by the irrep of the left
   % block. 
   % The Schmidt tensors have just one irrep label: the irrep of the left
   % block.
   % - LSMPS is fully left canonical by default.
   % - if obj.schmidt_list{end} is not just "1":  ==> infinite chain MPS
   %   (e.g. for iTEBD)

   obj.type = 'LSMPS';
   obj.matrices = cell(1,chain_length);
   obj.schmidt_list = cell(1,(chain_length));
   for pos = 1:chain_length
       obj.matrices{pos} = NAtensor({'t_left', 'tau', 'alpha', 't_right'},{'i','i','i','o'},{[1],[2],[1,2,3],[3]},no_of_symmetries);
       obj.schmidt_list{pos} = NAtensor({'t_left','t_right'},{'i','o'},{[1],[1]},no_of_symmetries);
   end
   obj.chain_length = chain_length;
   obj.no_of_symmetries = no_of_symmetries;
end

