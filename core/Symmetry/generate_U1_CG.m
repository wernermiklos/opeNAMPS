function CGtensor = generate_U1_CG(r1list, r2list, Rlist)
%GENERATE_U1_CG Generates the Clebsches of the U1 symmetry.
%   Detailed explanation goes here
    epsilon = 10^(-10);
    CGtensor = NAtensor({'m1', 'm2', 'M', 'alpha'}, {'i', 'i', 'o', 'o'} , ...
        {[1], [2], [3], [1, 2, 3]}, 1);
    
    for r1 = r1list
        for r2 = r2list
            % r1 and r2 are incoming representation labels, R is 'sum' of spins
            Q1 = (-1)^(r1 - 1)*floor(r1/2);
            Q2 = (-1)^(r2 - 1)*floor(r2/2);
            Q = Q1 + Q2;
            if abs(Q) < epsilon
                R = 1;
            elseif Q > 0
                R = 2*Q + 1;
            else
                R = - 2*Q;
            end
            if any(abs(Rlist-R) < epsilon)
                CGtensor = NTset_block(CGtensor,{{'m1',1},{'m2',1},{'M',1}}, {r1, r2, R}, ...
                        {'m1', 'm2', 'M', 'alpha'}, 1);
            end
        end
    end
end

