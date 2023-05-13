function out = LSMPS_set_schmidt(obj,pos, newschmidt)
            % Sets MPS matrix
            out = obj;
            out.schmidt_list{pos} = NAtensor({'t_left','t_right'},{'i','o'},{[1],[1]},obj.no_of_symmetries);
            out.schmidt_list{pos} = NTadd(out.schmidt{pos}, newschmidt);
end
