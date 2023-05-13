function out = LNTtoNT(obj)
            % Returns the tensor in full NAtensor format. (Be careful, it
            % may needs much-much memory!)
            out = obj.sub_tensors{1};
            for i = 2:obj.no_of_symmetries
                out = NTirrep_kron(out,obj.sub_tensors{i});
            end


end

