function CGtensor = generate_SU2_CG(r1list, r2list, Rlist,varargin)
%GENERATE_SU2 Generates the Clebsches of the SU2 symmetry.
%   Detailed explanation goes here
    epsilon = 10^(-10);
    parameters = struct();
    parameters.INFO = false;
    parameters = parameter_updater(parameters,varargin);
    info = parameters.INFO;

    CGtensor = NAtensor({'m1', 'm2', 'M', 'alpha'}, {'i', 'i', 'o', 'o'} , ...
        {[1], [2], [3], [1, 2, 3]}, 1);
    
    for r1 = r1list
        for r2 = r2list
            % r1 and r2 are incoming representation labels, R is 'sum' of spins
            j1 = (r1 - 1) / 2;
            j2 = (r2 - 1) / 2;
            ClebschBlocks = cell(1,r1+r2);   %Too many cells, but will work.
            for R = Rlist
                if info
                    disp([r1,r2,R])
                end
                J = (R - 1) / 2 ;
                if (J >= abs(j1-j2)) && (J <= j1 + j2) && (mod(r1+r2-1,2) == mod(R,2))
                    ClebschBlocks{R} = zeros(r1, r2, R, 1);
                    for m1 = (-j1):j1
                        i1 = j1 + 1 - m1;
                        for m2 = (-j2):j2
                            i2 = j2 + 1 - m2;
                            if (m1+m2 >= -J) && (m1+m2 <= J)
                                M = m1 + m2;
                                I = J + 1 - M;
                                ClebschBlocks{R}(i1,i2,I,1) = ClebschGordan_SU2(j1,j2,J,m1,m2,M);    %This is the function I stole from Ors.



                            end
                        end
                    end
                    CGtensor = NTset_block(CGtensor, {{'m1',1},{'m2',1},{'M',1}}, {r1, r2, R}, ...
                        {'m1', 'm2', 'M', 'alpha'}, ClebschBlocks{R});
                end
            end
        end
    end

                                           

    
end

