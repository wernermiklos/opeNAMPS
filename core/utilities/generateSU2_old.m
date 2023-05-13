function [IrrepDimensions,FusionMatrix,ClebschTensor] = generateSU2(MaxIrrepIndex)
%GENERATESU2 Summary of this function goes here
%   Detailed explanation goes here
FusionMatrix = zeros(MaxIrrepIndex,MaxIrrepIndex,MaxIrrepIndex);
for irrep1 = 1:MaxIrrepIndex
    for irrep2 = 1:MaxIrrepIndex
        for j = (abs(irrep2-irrep1) + 1):2:(min((irrep1+irrep2),MaxIrrepIndex))
            FusionMatrix(irrep1,irrep2,j) = 1;
        end
    end
end

IrrepDimensions = 1:MaxIrrepIndex;
irreprule = @(x) IrrepDimensions(x);
fusionrule = @(in) FusionMatrix(in(1),in(2),in(3));

ClebschTensor = NAtensor({'j1','j2','J'},{'m1','m2','M','eta'},{{'j1'},{'j2'},{'J'},{'j1','j2','J'}},MaxIrrepIndex);
ClebschTensor = activate_all(ClebschTensor,{'m1','m2','M','eta'},{irreprule,irreprule,irreprule,fusionrule});

for j1 = 0:0.5:((MaxIrrepIndex - 1)/2.0)
    irrep1 = 2*j1 + 1;
    for j2 = 0:0.5:((MaxIrrepIndex - 1)/2.0)
        irrep2 = 2*j2 + 1;
        ClebschBlocks = cell(2*(j1+j2)+1,1);
        for J = (j1+j2):(-1):(abs(j2-j1))
            IRREP = 2*J + 1;
            ClebschBlocks{IRREP} = zeros(irrep1,irrep2,IRREP);
            if J == j1+j2
                ClebschBlocks{IRREP}(1,1,1) = 1;
            else
                m1min = max(J-j2,-j1);
                m1max = min(J+j2,j1);
                A = zeros(j1+j2-J+1,j1+j2-J);
                for JTilde = (J+1):1:(j1+j2)
                    IrrepTilde = 2*JTilde + 1;
                    column = JTilde - J;
                    for m1=m1min:1:m1max
                        row = m1-m1min+1;
                        indexm1 =  j1 + 1 - m1;
                        indexm2 =  j2 + 1 - (J-m1);
                        indexM = JTilde-J+1;
                        A(row,column) = ClebschBlocks{IrrepTilde}(indexm1,indexm2,indexM);
                    end
                end
                b = findOrthogonal(A);
                for m1=m1min:1:m1max
                    row = m1-m1min+1;
                    indexm1 =  j1 + 1 - m1;
                    indexm2 =  j2 + 1 - (J-m1);
                    indexM = 1;
                    ClebschBlocks{IRREP}(indexm1,indexm2,indexM) = b(row);
                end
            end
            for M = (J-1):(-1):(-J)
                indexM = J + 1 - M;
                for m1 = (-j1):1:j1
                    indexm1 = j1 + 1 - m1;
                    for m2 = (-j2):1:j2
                        indexm2 = j2 + 1 - m2;
                        if m1 ~= j1
                            ClebschBlocks{IRREP}(indexm1,indexm2,indexM) =...
                                 ClebschBlocks{IRREP}(indexm1,indexm2,indexM) + ...
                                 ClebschBlocks{IRREP}(indexm1 - 1,indexm2,indexM - 1) * ...
                                 sqrt(j1*(j1+1)-m1*(m1+1))/sqrt(J*(J+1)-M*(M+1));
                        end
                        if m2 ~= j2
                            ClebschBlocks{IRREP}(indexm1,indexm2,indexM) =...
                                 ClebschBlocks{IRREP}(indexm1,indexm2,indexM) + ...
                                 ClebschBlocks{IRREP}(indexm1,indexm2 - 1,indexM - 1) * ...
                                 sqrt(j2*(j2+1)-m2*(m2+1))/sqrt(J*(J+1)-M*(M+1));
                        end
                    end
                end
            end
            if IRREP <= MaxIrrepIndex
                ClebschTensor = set_data_block(ClebschTensor,{'j1','j2','J'},[irrep1,irrep2,IRREP],...
                    {'m1','m2','M','eta'},ClebschBlocks{IRREP});
            end
        end
    end
end
                        
                        
            
            




end

