function out = SITE_generate_SU2xU1_spinhalf_fermion_open(SU2sym,U1sym)
%SITE_GENERATE_S1_SITE_SU2 creates a spin half fermionic site
%   Detailed explanation goes here


%spin-charge, pasku cikke alapjan
first_rep = [1,4];
second_rep = [1,1];
third_rep = [1,5];
fourth_rep = [2,2];
fifth_rep = [2,3];
sixth_rep = [3,1];
fdag_rep = [2,3];
f_rep = [2,2];
ftilddag_rep = [2,2];
ftild_rep = [2,3];
triv_rep = [1,1];

out = SITE_create({SU2sym,U1sym},{{first_rep,1},{second_rep,3},{third_rep,1},{fourth_rep,4},{fifth_rep,4},{sixth_rep,3}},'FERMIONS',true);

%rep: hova honnan
cdag_u_mx41=zeros(2,2,1,1);
cdag_u_mx41(2,1,1,1)=1;
cdag_u_mx52=zeros(2,2,3,1);
cdag_u_mx52(1,1,1,1)=1;
cdag_u_mx52(2,1,2,1)=-1/sqrt(2);
cdag_u_mx24=zeros(3,1,2,2);
cdag_u_mx24(2,1,1,1)=1/sqrt(2);
cdag_u_mx24(3,1,2,2)=1;
cdag_u_mx64=zeros(1,3,2,2);
cdag_u_mx64(1,2,1,1)=1/sqrt(2);
cdag_u_mx64(1,3,1,2)=1;
cdag_u_mx35=zeros(1,1,2,2);
cdag_u_mx35(1,1,1,2)=1;
cdag_u_mx56=zeros(2,2,1,3);
cdag_u_mx56(2,2,1,1)=1;
cdag_u_mx56(2,1,1,2)=1/sqrt(2);
cdag_u = NAtensor({'tau','mu','tau~','mu~'},{'o','o','i','i'},{[1],[1],[2],[2]},2);
cdag_u = NTset_block(cdag_u,{{'tau~',1},{'tau',1}},{fourth_rep,first_rep},{'tau~','mu~','tau','mu'},cdag_u_mx41);
cdag_u = NTset_block(cdag_u,{{'tau~',1},{'tau',1}},{fifth_rep,second_rep},{'tau~','mu~','tau','mu'},cdag_u_mx52);
cdag_u = NTset_block(cdag_u,{{'tau~',1},{'tau',1}},{second_rep,fourth_rep},{'tau~','mu~','tau','mu'},cdag_u_mx24);
cdag_u = NTset_block(cdag_u,{{'tau~',1},{'tau',1}},{sixth_rep,fourth_rep},{'tau~','mu~','tau','mu'},cdag_u_mx64);
cdag_u = NTset_block(cdag_u,{{'tau~',1},{'tau',1}},{third_rep,fifth_rep},{'tau~','mu~','tau','mu'},cdag_u_mx35);
cdag_u = NTset_block(cdag_u,{{'tau~',1},{'tau',1}},{fifth_rep,sixth_rep},{'tau~','mu~','tau','mu'},cdag_u_mx56);


cdag_d_mx41=zeros(2,2,1,1);
cdag_d_mx41(2,2,1,1)=1;
cdag_d_mx52=zeros(2,2,3,1);
cdag_d_mx52(1,2,1,1)=1;
cdag_d_mx52(2,2,2,1)=-1/sqrt(2);
cdag_d_mx24=zeros(3,1,2,2);
cdag_d_mx24(2,1,1,2)=-1/sqrt(2);
cdag_d_mx24(3,1,2,1)=-1;
cdag_d_mx64=zeros(1,3,2,2);
cdag_d_mx64(1,1,1,1)=1;
cdag_d_mx64(1,2,1,2)=1/sqrt(2);
cdag_d_mx35=zeros(1,1,2,2);
cdag_d_mx35(1,1,1,1)=-1;
cdag_d_mx56=zeros(2,2,1,3);
cdag_d_mx56(2,1,1,3)=-1;
cdag_d_mx56(2,2,1,2)=-1/sqrt(2);
cdag_d = NAtensor({'tau','mu','tau~','mu~'},{'o','o','i','i'},{[1],[1],[2],[2]},2);
cdag_d = NTset_block(cdag_d,{{'tau~',1},{'tau',1}},{fourth_rep,first_rep},{'tau~','mu~','tau','mu'},cdag_d_mx41);
cdag_d = NTset_block(cdag_d,{{'tau~',1},{'tau',1}},{fifth_rep,second_rep},{'tau~','mu~','tau','mu'},cdag_d_mx52);
cdag_d = NTset_block(cdag_d,{{'tau~',1},{'tau',1}},{second_rep,fourth_rep},{'tau~','mu~','tau','mu'},cdag_d_mx24);
cdag_d = NTset_block(cdag_d,{{'tau~',1},{'tau',1}},{sixth_rep,fourth_rep},{'tau~','mu~','tau','mu'},cdag_d_mx64);
cdag_d = NTset_block(cdag_d,{{'tau~',1},{'tau',1}},{third_rep,fifth_rep},{'tau~','mu~','tau','mu'},cdag_d_mx35);
cdag_d = NTset_block(cdag_d,{{'tau~',1},{'tau',1}},{fifth_rep,sixth_rep},{'tau~','mu~','tau','mu'},cdag_d_mx56);



ctilddag_u_mx42=zeros(2,2,3,1);
ctilddag_u_mx42(1,1,1,1)=1;
ctilddag_u_mx42(2,2,2,1)=1/sqrt(2);
ctilddag_u_mx53=zeros(2,2,1,1);
ctilddag_u_mx53(2,2,1,1)=1;
ctilddag_u_mx14=zeros(1,1,2,2);
ctilddag_u_mx14(1,1,1,2)=1;
ctilddag_u_mx25=zeros(3,1,2,2);
ctilddag_u_mx25(2,1,1,1)=-1/sqrt(2);
ctilddag_u_mx25(3,1,2,1)=1;
ctilddag_u_mx65=zeros(1,3,2,2);
ctilddag_u_mx65(1,1,1,2)=-1;
ctilddag_u_mx65(1,2,1,1)=-1/sqrt(2);
ctilddag_u_mx46=zeros(2,2,1,3);
ctilddag_u_mx46(2,1,1,3)=-1;
ctilddag_u_mx46(2,2,1,2)=-1/sqrt(2);
ctilddag_u = NAtensor({'tau','mu','tau~','mu~'},{'o','o','i','i'},{[1],[1],[2],[2]},2);
ctilddag_u = NTset_block(ctilddag_u,{{'tau~',1},{'tau',1}},{fourth_rep,second_rep},{'tau~','mu~','tau','mu'},ctilddag_u_mx42);
ctilddag_u = NTset_block(ctilddag_u,{{'tau~',1},{'tau',1}},{fifth_rep,third_rep},{'tau~','mu~','tau','mu'},ctilddag_u_mx53);
ctilddag_u = NTset_block(ctilddag_u,{{'tau~',1},{'tau',1}},{first_rep,fourth_rep},{'tau~','mu~','tau','mu'},ctilddag_u_mx14);
ctilddag_u = NTset_block(ctilddag_u,{{'tau~',1},{'tau',1}},{second_rep,fifth_rep},{'tau~','mu~','tau','mu'},ctilddag_u_mx25);
ctilddag_u = NTset_block(ctilddag_u,{{'tau~',1},{'tau',1}},{sixth_rep,fifth_rep},{'tau~','mu~','tau','mu'},ctilddag_u_mx65);
ctilddag_u = NTset_block(ctilddag_u,{{'tau~',1},{'tau',1}},{fourth_rep,sixth_rep},{'tau~','mu~','tau','mu'},ctilddag_u_mx46);


ctilddag_d_mx42=zeros(2,2,3,1);
ctilddag_d_mx42(1,2,1,1)=1;
ctilddag_d_mx42(2,1,2,1)=1/sqrt(2);
ctilddag_d_mx53=zeros(2,2,1,1);
ctilddag_d_mx53(2,1,1,1)=1;
ctilddag_d_mx14=zeros(1,1,2,2);
ctilddag_d_mx14(1,1,1,1)=-1;
ctilddag_d_mx25=zeros(3,1,2,2);
ctilddag_d_mx25(2,1,1,2)=1/sqrt(2);
ctilddag_d_mx25(3,1,2,2)=-1;
ctilddag_d_mx65=zeros(1,3,2,2);
ctilddag_d_mx65(1,3,1,1)=-1;
ctilddag_d_mx65(1,2,1,2)=-1/sqrt(2);
ctilddag_d_mx46=zeros(2,2,1,3);
ctilddag_d_mx46(2,1,1,2)=1/sqrt(2);
ctilddag_d_mx46(2,2,1,1)=1;
ctilddag_d = NAtensor({'tau','mu','tau~','mu~'},{'o','o','i','i'},{[1],[1],[2],[2]},2);
ctilddag_d = NTset_block(ctilddag_d,{{'tau~',1},{'tau',1}},{fourth_rep,second_rep},{'tau~','mu~','tau','mu'},ctilddag_d_mx42);
ctilddag_d = NTset_block(ctilddag_d,{{'tau~',1},{'tau',1}},{fifth_rep,third_rep},{'tau~','mu~','tau','mu'},ctilddag_d_mx53);
ctilddag_d = NTset_block(ctilddag_d,{{'tau~',1},{'tau',1}},{first_rep,fourth_rep},{'tau~','mu~','tau','mu'},ctilddag_d_mx14);
ctilddag_d = NTset_block(ctilddag_d,{{'tau~',1},{'tau',1}},{second_rep,fifth_rep},{'tau~','mu~','tau','mu'},ctilddag_d_mx25);
ctilddag_d = NTset_block(ctilddag_d,{{'tau~',1},{'tau',1}},{sixth_rep,fifth_rep},{'tau~','mu~','tau','mu'},ctilddag_d_mx65);
ctilddag_d = NTset_block(ctilddag_d,{{'tau~',1},{'tau',1}},{fourth_rep,sixth_rep},{'tau~','mu~','tau','mu'},ctilddag_d_mx46);

c_u=NTrename_legs(NTconj(cdag_u),{'tau','tau~','mu','mu~'},{'tau~','tau','mu~','mu'});
c_d=NTrename_legs(NTconj(cdag_d),{'tau','tau~','mu','mu~'},{'tau~','tau','mu~','mu'});
ctild_u=NTrename_legs(NTconj(ctilddag_u),{'tau','tau~','mu','mu~'},{'tau~','tau','mu~','mu'});
ctild_d=NTrename_legs(NTconj(ctilddag_d),{'tau','tau~','mu','mu~'},{'tau~','tau','mu~','mu'});


ph_mx11=zeros(1,1,1,1);
ph_mx11(1,1,1,1)=1;
ph_mx22=zeros(3,1,3,1);
ph_mx22(1,1,1,1)=1;
ph_mx22(2,1,2,1)=1;
ph_mx22(3,1,3,1)=1;
ph_mx33=zeros(1,1,1,1);
ph_mx33(1,1,1,1)=1;
ph_mx44=zeros(2,2,2,2);
ph_mx44(1,1,1,1)=-1;
ph_mx44(1,2,1,2)=-1;
ph_mx44(2,1,2,1)=-1;
ph_mx44(2,2,2,2)=-1;
ph_mx55=zeros(2,2,2,2);
ph_mx55(1,1,1,1)=-1;
ph_mx55(1,2,1,2)=-1;
ph_mx55(2,1,2,1)=-1;
ph_mx55(2,2,2,2)=-1;
ph_mx66=zeros(1,3,1,3);
ph_mx66(1,1,1,1)=1;
ph_mx66(1,2,1,2)=1;
ph_mx66(1,3,1,3)=1;
ph = NAtensor({'tau','mu','tau~','mu~'},{'o','o','i','i'},{[1],[1],[2],[2]},2);
ph = NTset_block(ph,{{'tau~',1},{'tau',1}},{first_rep,first_rep},{'tau~','mu~','tau','mu'},ph_mx11);
ph = NTset_block(ph,{{'tau~',1},{'tau',1}},{second_rep,second_rep},{'tau~','mu~','tau','mu'},ph_mx22);
ph = NTset_block(ph,{{'tau~',1},{'tau',1}},{third_rep,third_rep},{'tau~','mu~','tau','mu'},ph_mx33);
ph = NTset_block(ph,{{'tau~',1},{'tau',1}},{fourth_rep,fourth_rep},{'tau~','mu~','tau','mu'},ph_mx44);
ph = NTset_block(ph,{{'tau~',1},{'tau',1}},{fifth_rep,fifth_rep},{'tau~','mu~','tau','mu'},ph_mx55);
ph = NTset_block(ph,{{'tau~',1},{'tau',1}},{sixth_rep,sixth_rep},{'tau~','mu~','tau','mu'},ph_mx66);

n_u=NTdot(cdag_u, c_u, {'tau','mu'},{'tau~','mu~'});
n_d=NTdot(cdag_d, c_d, {'tau','mu'},{'tau~','mu~'});
ntild_u=NTdot(ctilddag_u, ctild_u, {'tau','mu'},{'tau~','mu~'});
ntild_d=NTdot(ctilddag_d, ctild_d, {'tau','mu'},{'tau~','mu~'});

n = NTadd(n_u, n_d);
ntild = NTadd(ntild_u, ntild_d);

c_u__ctild_u = NTdot(c_u, ctild_u, {'tau','mu'},{'tau~','mu~'});
c_d__ctild_d = NTdot(c_d, ctild_d, {'tau','mu'},{'tau~','mu~'});


%Hubbard coupling: fdag_u fdag_d f_d f_u
hub = NTdot(n_u,n_d,{'tau','mu'}, {'tau~','mu~'});
hubtild = NTdot(ntild_u,ntild_d,{'tau','mu'}, {'tau~','mu~'});

c__ctild=NTadd(c_u__ctild_u,c_d__ctild_d);


out = SITE_define_tensor_operator(out,'fdag',{SU2sym,U1sym},fdag_rep,{cdag_u,cdag_d},'FERMIONICITY',-1);
out = SITE_define_tensor_operator(out,'f',{SU2sym,U1sym},f_rep,{c_d,NTneg(c_u)},'FERMIONICITY',-1);
out = SITE_define_tensor_operator(out,'ftilddag',{SU2sym,U1sym},ftilddag_rep,{ctilddag_d,ctilddag_u},'FERMIONICITY',-1);
out = SITE_define_tensor_operator(out,'ftild',{SU2sym,U1sym},ftild_rep,{ctild_u,NTneg(ctild_d)},'FERMIONICITY',-1); 
out = SITE_define_tensor_operator(out,'ph',{SU2sym,U1sym},triv_rep,{ph},'FERMIONICITY',+1);
out = SITE_define_tensor_operator(out,'n',{SU2sym,U1sym},triv_rep,{n},'FERMIONICITY',+1);
out = SITE_define_tensor_operator(out,'ntild',{SU2sym,U1sym},triv_rep,{ntild},'FERMIONICITY',+1);
out = SITE_define_tensor_operator(out,'hub',{SU2sym,U1sym},triv_rep,{hub},'FERMIONICITY',+1);
out = SITE_define_tensor_operator(out,'hubtild',{SU2sym,U1sym},triv_rep,{hubtild},'FERMIONICITY',+1);
out = SITE_define_tensor_operator(out,'c__ctild',{SU2sym,U1sym},triv_rep,{c__ctild},'FERMIONICITY',+1);


end
