function out = LNTunprime_all_names(obj)
            out = LNTcopy(obj);
            for symID = 1:obj.no_of_symmetries
                out.sub_tensors{symID} = NTunprime_all_names(out.sub_tensors{symID});
            end
            out.leg_names = out.sub_tensors{1}.leg_names;

end

