function out = NTsubtr(A, B)
% Subtracts two NAtensors ("obj" and "other"),  if their structure 
           % (leg_names, leg_types, dependency structure, block shapes) are consistent.
           %
           % It is also allowed to add a tensor to the 0 number, then the
           % nonzero tensor is just copied.
           % ----
           % obj:    - NAtensor
           % other:  - the other NAtensor
           % ---
           % out = obj - other.
    out = NTadd(A, NTneg(B));
end

