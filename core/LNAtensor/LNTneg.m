function out = LNTneg(obj)
           % Overloading "-a"
           if strcmp(obj.type,'LNAtensor')
               out = LNTmult(obj,(-1.0));
           elseif strcmp(obj.type,'NAtensor')
               out = NTmult(obj,(-1.0));
           else
               error('LNTneg works for LNAtensors and NAtensors');
           end
end

