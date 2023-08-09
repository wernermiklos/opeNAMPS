NO_OF_SYMS = 2;

% Definition of an NA-tensor
T = NAtensor({'a','b','c','eta'}, ...       %leg names
             {'i','i','o','o'}, ...         %leg types
             {[1],[2],[3],[1,2,3]}, ...     %dependencies
             NO_OF_SYMS);                   %number of symmetries

% Setting a block
T = NTset_block(T, ...
                {{'a',1},{'b',1}, {'c',1}}, ...
                {[2,3] ,[3,4], [4,5]}, ...
                {'a','b','c','eta'}, ...
                ones(4,6,8,2));

% alternative, but equivalent block setting
%T = NTset_block(T, ...
%                {{'b',1},{'eta',3}, {'eta',1}}, ...
%                {[3,4] ,[4,5], [2,3]}, ...
%                {'eta','b','a','c'}, ...
%                ones(2,6,4,8));

twoT = NTmult(T,2);
TplusT = NTadd(T,T);
TminustwoT = NTsubtr(T,twoT);
minusT = NTneg(T);
conjT = NTconj(T);

%NTprime_all_legs: adds '~' to leg names
conjT_primed = NTprime_all_names(conjT);
TconjT = NTdot(T,conjT_primed,{'a'},{'a~'});


