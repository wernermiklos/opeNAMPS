function out = cartprod_fromcell(incell)
%CARTPROD_FROMCELL Works for row vectors
    if length(incell) == 1
        out = incell{1};
        return
    end
    tmp = [kron(incell{1},ones(1,size(incell{2},2))); ...
           kron(ones(1,size(incell{1},2)),incell{2})];
    for i = 3:length(incell)
        tmp = [kron(tmp,ones(1,size(incell{i},2))); ...
               kron(ones(1,size(tmp,2)),incell{i})];
    end
    out = tmp;

end

