function out = SITE_generate_SU2_spin_site(SU2sym, S)
%SITE_GENERATE_SU2_SPIN_SITE generates a spin site (spin length: S)
%   Detailed explanation goes here
 site_rep = 2*S + 1;
 spinop_rep = 3;   % spin is a vector operator
 triv_rep = 1;


 out = SITE_create({SU2sym},{{site_rep,1}});

 SpMatrix = zeros(2*S+1,2*S+1);
 SzMatrix = zeros(2*S+1,2*S+1);
 for i = 1:(2*S+1)
     M = S - i + 1;
     SzMatrix(i,i) = M;
     if i ~= 1
        SpMatrix(i-1,i) = sqrt(S*(S+1) - M*(M+1));
     end
 end
 SmMatrix = SpMatrix';

 Sp = NAtensor({'tau','mu','tau~','mu~'},{'o','o','i','i'},{[1],[1],[2],[2]},1);   % tau legs are dummy one dim legs in spin sites.
 Sp = NTset_block(Sp,{{'tau~',1},{'tau',1}},{site_rep,site_rep},{'mu~','mu','tau~','tau'},SpMatrix);
  
 Sz = NAtensor({'tau','mu','tau~','mu~'},{'o','o','i','i'},{[1],[1],[2],[2]},1);   % tau legs are dummy one dim legs in spin sites.
 Sz = NTset_block(Sz,{{'tau~',1},{'tau',1}},{site_rep,site_rep},{'mu~','mu','tau~','tau'},SzMatrix);

 Sm = NAtensor({'tau','mu','tau~','mu~'},{'o','o','i','i'},{[1],[1],[2],[2]},1);   % tau legs are dummy one dim legs in spin sites.
 Sm = NTset_block(Sm,{{'tau~',1},{'tau',1}},{site_rep,site_rep},{'mu~','mu','tau~','tau'},SmMatrix);

 out = SITE_define_tensor_operator(out,'S',{SU2sym},spinop_rep, ...
                                   {NTmult(Sp,-1/sqrt(2)), ...
                                    Sz, ...
                                    NTmult(Sm,1/sqrt(2))});

end

