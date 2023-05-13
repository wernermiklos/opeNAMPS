function out = LNTirrep_kron(a, b, varargin)
            % Returns irrep_kron of a and b.
            if strcmp(a.type,'LNAtensor')
                obj = a;
                other = b;
                reversed = false;
            elseif strcmp(b.type,'LNAtensor')
                obj = b;
                other = a;
                reversed = true;
            elseif strcmp(a.type,'NAtensor') && strcmp(b.type,'NAtensor')
                out = NTirrep_kron(a,b,'LNAtensor');
                return;
            else
                error('a and b must be NAtensor or LNAtensor')
            end
            if ~strcmp(other.type,'NAtensor') && ~strcmp(other.type,'LNAtensor')
                error('Error: other must be an NAtensor too.')
            end
            
            if strcmp(other.type,'NAtensor')
                if other.no_of_symmetries ~= 1
                    error('If other is NAtensor it must have no_of_symmetries = 1')
                end
                out = LNAtensor(obj.leg_names,obj.leg_types,obj.dependencies,obj.no_of_symmetries + 1);
                if reversed
                    objsymIDrange = 2:(out.no_of_symmetries);
                    othersymID = 1;
                else
                    objsymIDrange = 1:(out.no_of_symmetries-1);
                    othersymID = out.no_of_symmetries;
                end
                out = LNTset_subtensor(out,othersymID,other);
                for symID = 1:obj.no_of_symmetries
                    out = LNTset_subtensor(out, objsymIDrange(symID),obj.sub_tensors{symID});
                end
            else
                out = LNAtensor(obj.leg_names,obj.leg_types,obj.dependencies,obj.no_of_symmetries + other.no_of_symmetries);
                for symID = 1:obj.no_of_symmetries
                    out = LNTset_subtensor(out,symID,obj.sub_tensors{symID});
                end
                for othersymID = 1:other.no_of_symmetries
                    symID = obj.no_of_symmetries + othersymID;
                    out = LNTset_subtensor(out,symID,other.sub_tensors{othersymID});
                end
            end


end

