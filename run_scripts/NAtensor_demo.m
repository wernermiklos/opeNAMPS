NO_OF_SYMS = 2;

% Definition of an NA-tensor (first example: tensor with the structure of a
%                                            simple abelian MPS matrix)
M1 = NAtensor({'t_in','tau','t_out'}, ...    %leg names
              {'i','i','o'}, ...             %leg types
              {[1],[2],[3]}, ...             %dependencies (in this simple case each leg has its quantum number)
              NO_OF_SYMS);                   %number of symmetries (length of 'Gamma' indices)
M2 = NAtensor({'t_in','tau','t_out'}, ...    %leg names
              {'i','i','o'}, ...             %leg types
              {[1],[2],[3]}, ...             %dependencies (in this simple case each leg has its quantum number)
              NO_OF_SYMS);                   %number of symmetries (length of 'Gamma' indices)

% Setting a block
M1 = NTset_block(M1, ...
                 {{'t_in',1},{'tau',1}, {'t_out',1}}, ...   % order of irrep labels
                 {[2,3] ,[3,4], [4,5]}, ...                 % value of irrep labels
                 {'t_in','tau','t_out'}, ...                % order of legs
                 rand(2,3,4));                              % some random tensorblock.
M1 = NTset_block(M1, ...
                 {{'t_in',1},{'tau',1}, {'t_out',1}}, ...
                 {[1,4] ,[2,3], [4,5]}, ...
                 {'t_in','tau','t_out'}, ...
                 rand(1,5,4));        % careful, the irrep of t_out is the same as in the block before: 
                                      % the dimension of t_out must also be the same!

M2 = NTset_block(M2, ...
                 {{'t_in',1},{'tau',1}, {'t_out',1}}, ...   % order of irrep labels
                 {[4,5] ,[2,3], [5,6]}, ...                 % value of irrep labels
                 {'t_in','tau','t_out'}, ...                % order of legs
                 rand(4,5,5));                              % some random tensorblock.
M2 = NTset_block(M2, ...
                 {{'t_in',1},{'tau',1}, {'t_out',1}}, ...
                 {[4,5] ,[3,4], [2,3]}, ...
                 {'t_in','tau','t_out'}, ...
                 rand(4,3,2));        % careful, the irrep of t_in is the same as in the block before: 
                                      % the dimension of t_in must also be the same!

% We now want to contract the 't_out' leg of M1 with the 't_in' leg of M2.
% However, both tensors have uncontracted legs with the name 'tau', i.e. we
% have to resolve the name-redundancy.
% Three solutions are presented:
% i.) rename the legs by hand
M1newnames = NTrename_legs(M1,{'tau'},{'tau1'});
M2newnames = NTrename_legs(M2,{'tau'},{'tau2'});
M1dotM2_renamed_by_hand = NTdot(M1newnames,M2newnames,{'t_out'},{'t_in'});
% ii.) prime the leg names of M2;
M2primednames = NTprime_all_names(M2);
M1dotM2_primedM2 = NTdot(M1,M2primednames,{'t_out'},{'t_in~'});
% iii.) perform renaming on-the-fly by adding a 5th argument to NTdot
M1dotM2_onthefly = NTdot(M1,M2,{'t_out'},{'t_in'},{{},{{'tau','tau2'}}});   % we only rename 'tau' leg of M2.


%Simple algebraic manipulations

twoM1 = NTmult(M1,2);   % 2*M1
M1plusM2 = NTadd(M1,M2);  % M1 + M2 
M2minus_twoM1_1 = NTsubtr(M2,twoM1);  % M2 - 2*M1 (first way)
M2minus_twoM1_2 = NTadd(M2,NTmult(M1,-2));   % M2 - 2*M1 (second way)
minusM2 = NTneg(M2);    % -M2
conjM1 = NTconj(M1);    % conjugation reverts the direction of legs too.

% Definition of an NAtensor (second example: MPS matrix of a left canonical non-Abelian MPS)
A = NAtensor({'t_left','tau','alpha','t_right'}, ...    %leg names
              {'i','i','i','o'}, ...                    %leg types
              {[1],[2],[1,2,3],[3]}, ...                %dependencies (mind the triple dep of 'alpha', and the overlap between deps.)                                             
              NO_OF_SYMS);

% We set now two blocks, but in two different ways. Both are correct.
A = NTset_block(A,...
                {{'t_left',1},{'tau',1}, {'t_right',1}}, ...   % order of irrep labels
                {[2,3] ,[3,4], [4,5]}, ...                     % value of irrep labels
                {'alpha','t_left','tau','t_right'}, ...            % order of legs
                rand(1,2,3,4));                                % some random tensorblock.)
A = NTset_block(A,...
                {{'t_left',1},{'alpha',3}, {'alpha',2}}, ...   % order of irrep labels
                {[1,4] ,[4,5], [2,3]}, ...                     % value of irrep labels
                {'t_left','tau','alpha','t_right'}, ...            % order of legs
                rand(1,5,1,4));                                % some random tensorblock.)

% Demo of NTdot rules
conjAprimed = NTconj(NTprime_all_names(A));

AA1 = NTdot(A,conjAprimed,{'t_left'},{'t_left~'});
AA2 = NTdot(A,conjAprimed,{'tau','alpha'},{'tau~','alpha~'});
AA3 = NTdot(A,conjAprimed,{'t_left','tau','alpha'},{'t_left~','tau~','alpha~'});

