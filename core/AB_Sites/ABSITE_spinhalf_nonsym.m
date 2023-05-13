function out = ABSITE_spinhalf_nonsym()

%
%   Constructs a general site based on the Z1 symmetry (no symmetry...)
%
%
% 


Sz = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},1);
Sz = NTset_block(Sz,{{'tau',1},{'tau~',1}},{[1],[1]},{'tau~','tau'},[0.5 0; ...
                                                                     0  -0.5]);

Sx = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},1);
Sx = NTset_block(Sx,{{'tau',1},{'tau~',1}},{[1],[1]},{'tau~','tau'},[0 0.5; ...
                                                                     0.5 0]);

Sy = NAtensor({'tau','tau~'},{'o','i'},{[1],[2]},1);
Sy = NTset_block(Sy,{{'tau',1},{'tau~',1}},{[1],[1]},{'tau~','tau'},[0 -0.5i; ...
                                                                     0.5i 0]);

Sp = NTadd(Sx, NTmult(Sy,1i));
Sm = NTadd(Sx, NTmult(Sy,-1i));

out = ABSITE_create(1,{{[1],2}});

out.name = 'ABSITE_spinhalf_nonsym';

out = ABSITE_define_operator(out,'Sz',Sz,1);
out = ABSITE_define_operator(out,'S+',Sp,1);
out = ABSITE_define_operator(out,'S-',Sm,1);
out = ABSITE_define_operator(out,'Sx',Sx,1);
out = ABSITE_define_operator(out,'Sy',Sy,1);

end