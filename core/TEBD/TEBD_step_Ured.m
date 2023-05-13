function [Aout1, lout1, Aout2, truncerr, M_mult_actual] = TEBD_step_Ured( A1, A2, l2, U_red, M_mult,varargin)
%TEBD_STEP Summary of this function goes here
%   Detailed explanation goes here
    parameters = struct();
    parameters.EPS = 10^(-16);
    parameters.NORMALIZE = true;
    parameters = parameter_updater(parameters,varargin);

    eps = parameters.EPS;
    C = NTdot(A1,A2,{'t_right'},{'t_left'},...
            {{{'tau','tau_1'},{'alpha','alpha_1'}},...
             {{'tau','tau_2'},{'alpha','alpha_2'}}});
    newC = NTdot(C,U_red,...
               {'tau_1','alpha_1','tau_2','alpha_2'},...
               {'tau_1','alpha_1','tau_2','alpha_2'},...
               {{},{{'tau_1~','tau_1'},{'tau_2~','tau_2'},{'alpha_1~','alpha_1'},{'alpha_2~','alpha_2'}}});
    newTheta = NTdot(newC,l2,{'t_right'},{'t_left'});
    normsq_0 = NTdot(newTheta,NTconj(newTheta),...
                     {'t_left','tau_1','alpha_1','tau_2','alpha_2','t_right'}, ...
                     {'t_left','tau_1','alpha_1','tau_2','alpha_2','t_right'});
    [A1tmp, l1tmp, ~, normsq, M_mult_actual] = NTsvd(newTheta,...
                                                     {'t_left','tau_1','alpha_1'},...
                                                     {'tau_2','alpha_2','t_right'},...
                                                     {'t_right','t_left'},{'i','o'},M_mult, eps);
    l1tmp = NTrename_legs(l1tmp, {'t_right','t_left'},{'t_left','t_right'});
    A2tmp = NTdot(newC,NTconj(A1tmp),{'t_left','tau_1','alpha_1'},{'t_left','tau_1','alpha_1'},{{},{{'t_right','t_left'}}});
    A1tmp = NTrename_legs(A1tmp, {'tau_1','alpha_1'},{'tau','alpha'});
    A2tmp = NTrename_legs(A2tmp, {'tau_2','alpha_2'},{'tau','alpha'});
    Aout1 = A1tmp;
    if parameters.NORMALIZE
        lout1 = NTmult(l1tmp, (1/sqrt(normsq)));
        Aout2 = NTmult(A2tmp, (1/sqrt(normsq)));
    else
        lout1 = NTmult(l1tmp, (sqrt(normsq_0)/sqrt(normsq)));    %fix the norm to the untruncated value normsq_0
        Aout2 = NTmult(A2tmp, (sqrt(normsq_0)/sqrt(normsq)));
    end
    truncerr = (normsq_0-normsq)/normsq_0;
    
end

