function obj = ABSITE_define_operator(obj_in, ...
                                           name, ...                  % char array, the name of the operator (like 'fdag')
                                           operator, ...                % List of operators (the list index corresponds to the "M_op" index of the tensor operator)
                                           fermionicity)                 % Tolerance value (10^(-12) by default). It is used to check if the op is indeed a tensor op.
             if nargin == 3
                 fermionicity = 1;
             elseif nargin > 4
                 error('function has 3 or 4  arguments')
             end
             obj = obj_in;
             obj.operators(name) = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},obj.no_of_symmetries);
             obj.operators(name) = NTadd(obj.operators(name), operator);
             obj.fermionicities(name) = fermionicity;
end

