function  out = LNTadd(a, b)
            % Overloading obj + other. Works only if the two tensors are
            % consistent. Be careful, the result is a full NAtensor.
            if isequal(a,0)
                out = b;
                return
            elseif isequal(b,0)
                out = a;
                return
            end
            
            if strcmp(a.type,'LNAtensor') && strcmp(b.type,'LNAtensor')
                out = NTadd(LNTtoNT(a),LNTtoNT(b));
            elseif strcmp(a.type,'LNAtensor') && strcmp(b.type,'NAtensor')
                out = NTadd(LNTtoNT(a),b);
            elseif strcmp(a.type,'NAtensor') && strcmp(b.type,'LNAtensor')
                out = NTadd(a,LNTtoNT(b));
            elseif strcmp(a.type,'NAtensor') && strcmp(b.type,'NAtensor')
                out = NTadd(a,b);
            else
                error('LNTadd can add LNAtensors or NAtensors');
            end
            
end

