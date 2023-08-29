% To create a symmetry, we can use the built in symmetries (U1, SU2 & Z2)

Z2sym = generate_Z2_Symmetry();  % even: REP = 1; odd: REP = 2

U1sym = generate_U1_Symmetry(1:11,1:11,1:21);    %arguments: allowed replist in1, in2 and out  
% U1 representation table:
% Q   | 0 | -1 | +1 | -2 | +2 | -3 | +3 | ... 
% ----|---|----|----|----|----|----|----|------------
% REP | 1 |  2 |  3 |  4 |  5 |  6 |  7 | ...

SU2sym = generate_SU2_Symmetry(1:11,1:11,1:21);  % spin rep index: REP = 2*S+1

% In case of U1 we can generate the symmetry without the Clebsches. Abelian
% routines don't need the CG tensor.
U1symDummy = generate_U1_Symmetry([],[],[]);


% We can create a any symmetry from its CG tensor. E.G. re-create SU(2)
% from its CG tensor:
CGSU2 = generate_SU2_CG(1:11,1:11,1:21);

SU2sym_fromCG = SYM_create(CGSU2,'SU2 generated from CG');

%Compare the two SU2's. The selection rule / fusion rule / conj_reps fields
%are different! The built in is faster (the default is replaced). 

