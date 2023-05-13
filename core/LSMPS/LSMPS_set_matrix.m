function out = LSMPS_set_matrix(obj,pos, newmatrix)
            % Sets MPS matrix
            out = obj;
            out.matrices{pos} = NAtensor({'t_left','tau','alpha','t_right'},{'i','i','i','o'},{[1],[2],[1,2,3],[3]},obj.no_of_symmetries);
            out.matrices{pos} = NTadd(out.matrices{pos}, newmatrix);
end

