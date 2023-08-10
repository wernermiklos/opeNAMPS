function a = compare_SU2_CG_Pascu_Ors(r1list, r2list, Rlist,info)
%GENERATE_SU2 Generates the Clebsches of the SU2 symmetry.
%   Detailed explanation goes here
    epsilon = 10^(-10);
    if nargin == 3
        info = false;
    end
    a = zeros(0,8);
    for r1 = r1list
        for r2 = r2list
            % r1 and r2 are incoming representation labels, R is 'sum' of spins
            j1 = (r1 - 1) / 2;
            j2 = (r2 - 1) / 2;
            for R = Rlist
                if info
                    disp([r1,r2,R])
                end
                J = (R - 1) / 2 ;
                if (J >= abs(j1-j2)) && (J <= j1 + j2) && (mod(r1+r2-1,2) == mod(R,2))
                    for m1 = (-j1):j1
                        i1 = j1 + 1 - m1;
                        for m2 = (-j2):j2
                            i2 = j2 + 1 - m2;
                            if (m1+m2 >= -J) && (m1+m2 <= J)
                                M = m1 + m2;
                                I = J + 1 - M;
                                diff = clebschgordan(j1,m1,j2,m2,J,M) - ClebschGordan_SU2(j1,j2,J,m1,m2,M);    %This is the function I got from Pascu.
                                if abs(diff) > epsilon
                                    disp(abs(diff))
                                    a(end+1,:) = [j1,m1,j2,m2,J,M,clebschgordan(j1,m1,j2,m2,J,M), ClebschGordan_SU2(j1,j2,J,m1,m2,M)];
                                end


                            end
                        end
                    end
                end
            end
        end
    end

                                           

    
end

