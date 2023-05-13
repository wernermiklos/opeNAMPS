function CGtensor = generate_SU2_CG(r1list, r2list, Rlist)
%GENERATE_SU2 Generates the Clebsches of the SU2 symmetry.
%   Detailed explanation goes here
    epsilon = 10^(-10);
    CGtensor = NAtensor({'m1', 'm2', 'M', 'alpha'}, {'i', 'i', 'o', 'o'} , ...
        {[1], [2], [3], [1, 2, 3]}, 1);
    
    for r1 = r1list
        for r2 = r2list
            % r1 and r2 are incoming representation labels, R is 'sum' of spins
            j1 = (r1 - 1) / 2;
            j2 = (r2 - 1) / 2;
            ClebschBlocks = cell(1,r1+r2);   %Too many cells, but will work.
            for R = (r1+r2-1):-2:(abs(r1-r2)+1)
                J = (R - 1) / 2 ;
                ClebschBlocks{R} = zeros(r1, r2, R, 1);
                if abs(J - (j1 + j2)) < epsilon 
                    ClebschBlocks{R}(1, 1, 1, 1) = 1.0;
                else
                    subspace_dim = (j1 + j2 - J) + 1;
                    M = J;
                    newvec = zeros(1,subspace_dim);
                    others = zeros(subspace_dim, subspace_dim - 1);
                    i = 1;
                    for ind1 = 1:r1
                        for ind2 = 1:r2
                            m1 = j1-ind1+1;
                            m2 = j2-ind2+1;
                            if abs(m1 + m2 - M) < epsilon
                                for colind = 1:(subspace_dim - 1)
                                    Jtmp = j1 + j2 - colind + 1;
                                    Rtmp = 2*Jtmp + 1;
                                    IND = Jtmp-M+1;
                                    others(i, colind) = ClebschBlocks{Rtmp}(...
                                        ind1, ind2, IND, 1);
                                end
                                i = i + 1;
                            end
                        end
                    end
                    i = 1;
                    for ind1 = 1:r1
                        for ind2 = 1:r2
                            m1 = j1 - ind1 + 1;
                            m2 = j2 - ind2 + 1;
                            if abs(m1 + m2  - M) < epsilon
                                minormx = others([1:(i-1),(i+1):end],:);
                                newvec(i) = det(minormx) * (-1)^i;
                                ClebschBlocks{R}(ind1, ind2, 1, 1) = newvec(i);
                                i = i + 1;
                            end
                        end
                    end
                end
                for IND = 2:R
                    M = J - IND + 1;
                    for ind1 = 1:r1
                        for ind2 = 1:r2
                            m1 = j1-ind1+1;
                            m2 = j2-ind2+1;
                            if abs(m1 + m2 - M) < epsilon
                                if ind1 ~= 1 && ind2 ~= 1
                                ClebschBlocks{R}(ind1, ind2, IND, 1) = 1.0 / ...
                                    sqrt(J * (J + 1.) - M * (M + 1.)) * ...
                                    (sqrt(j1 * (j1 + 1.) - m1 * (m1 + 1)) * ...
                                    ClebschBlocks{R}(ind1 - 1, ind2, IND - 1, 1) + ...
                                    sqrt(j2 * (j2 + 1.) - m2 * (m2 + 1)) * ...
                                    ClebschBlocks{R}(ind1, ind2 - 1, IND - 1, 1));
                                end
                                if ind1 == 1 && ind2 ~= 1
                                ClebschBlocks{R}(ind1, ind2, IND, 1) = 1.0 / ...
                                    sqrt(J * (J + 1.) - M * (M + 1.)) * ...
                                    (sqrt(j2 * (j2 + 1.) - m2 * (m2 + 1)) * ...
                                    ClebschBlocks{R}(ind1, ind2 - 1, IND - 1, 1));
                                end
                                if ind1 ~= 1 && ind2 == 1
                                ClebschBlocks{R}(ind1, ind2, IND, 1) = 1.0 / ...
                                    sqrt(J * (J + 1.) - M * (M + 1.)) * ...
                                    (sqrt(j1 * (j1 + 1.) - m1 * (m1 + 1)) * ...
                                    ClebschBlocks{R}(ind1 - 1, ind2, IND - 1, 1));
                                end
                                
                            end
                        end
                    end
                end
                if any(abs(Rlist-R) < epsilon)
                    CGtensor = NTset_block(CGtensor, {{'m1',1},{'m2',1},{'M',1}}, {r1, r2, R}, ...
                        {'m1', 'm2', 'M', 'alpha'}, ClebschBlocks{R});
                end
            end
        end
    end

                                           

    
end

