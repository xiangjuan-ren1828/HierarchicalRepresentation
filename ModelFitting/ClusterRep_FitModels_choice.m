function [lliEst, paramsEst, lli_all, fn, fnDir] = ClusterRep_FitModels_choice(expID, subID, stim, resp, states, choiceId, delta_t_trials, assDist_mat, fitWord, ...
                                        expMode, posXY, respYes, mouseTraj, objAng_col, objDtrNo_col, miniDmat, nFit, Midx, folder, refit, transMat_input, nij_input, iteWord, trialWord, priorWord)     

% write by rxj @ 08/10/2021 based on ClusterRep_FitModels.m
% =================================================
% write by rxj @ 05/05/2021 based on SR_FitModels.m
% ------------Description------------
% dataCleanID: preprocess procedure (full length vs. only hamiltonian trials; take all effect into consideration vs. regress-out)
% expID:       which experiment (implicit, explicit, implicit+random)
% Midx:        model assumptions
% =================================================
% write by rxj @ 07/24/2020

if isequal(expMode, 'mouse')
    if isequal(trialWord, 'all')
        if isequal(Midx, 'MaxEntropy') && isequal(priorWord, 'InformPrior')
            fnDir = [folder, '/ClusterRep_ChoiceModel_Results/', expID, '/', Midx, '-ModelFits-', fitWord, '/InformPrior/'];
        else
            fnDir = [folder, '/ClusterRep_ChoiceModel_Results/', expID, '/', Midx, '-ModelFits-', fitWord, '/updates2/'];
        end
    elseif isequal(trialWord, 'Hamiltonian') || isequal(trialWord, 'Hamiltonian-binary') % Hamiltonian-binary: only trials with 2 options, one target and another distractor
        fnDir = [folder, '/ClusterRep_ChoiceModel_Results/', expID, '-', trialWord, '/', Midx, '-ModelFits-', fitWord, '/'];
    end
    
elseif isequal(expMode, 'key')
    if isequal(trialWord, 'all')
        if isequal(Midx, 'MaxEntropy') && isequal(priorWord, 'InformPrior')
            fnDir = [folder, '/ClusterRep_ChoiceModel_Results/NormalDist-stdFixNorm/', expID, '/', Midx, '-ModelFits/', '/InformPrior/'];
        else
            fnDir = [folder, '/ClusterRep_ChoiceModel_Results/NormalDist-stdFixNorm/', expID, '/', Midx, '-ModelFits/'];
        end
    elseif isequal(trialWord, 'Hamiltonian')
        fnDir = [folder, '/ClusterRep_ChoiceModel_Results/NormalDist', '-', trialWord, '/', expID, '/', Midx, '-ModelFits/'];
    end
end
if ~exist(fnDir, 'dir')
    mkdir(fnDir);
end
if isequal(iteWord, 'naive')       %% learn the transition matrix from zero
    fn = [fnDir, expID, '-', subID, '-rep', num2str(nFit), '.mat'];
elseif isequal(iteWord, 'trained')
    fn = [fnDir, expID, '-', iteWord, '-', subID, '-rep', num2str(nFit), '.mat'];
end

if exist(fn, 'file') && refit == 0
    load(fn, 'lliEst', 'paramsEst', 'lli_all');
else
    if isequal(Midx, 'MaxEntropy')
        theModel = @maxEntropy_model_choice;
        if isequal(expMode, 'mouse')
            LB = [0, 0]; % beta (inverse temperature), r1 (inverse temperateure of the choice model)
            UB = [Inf, Inf];
            if isequal(priorWord, 'InformPrior')
                LB = [LB, 0]; % the initial value of the matrix nij
                UB = [UB, Inf];
            end
        elseif isequal(expMode, 'key')
            LB = [0, -Inf]; % beta (inverse temperature), L (decision criterion)
            UB = [Inf, Inf];
            if isequal(priorWord, 'InformPrior')
                LB = [LB, 0]; % the initial value of the matrix nij
                UB = [UB, Inf];
            end
        end
        
    elseif isequal(Midx, 'SR')
        theModel = @SR_model_choice;
        if isequal(expMode, 'mouse')
            LB = [0, 0, 0]; % gamma (discount factor), alpha (learning rate), beta (inverse temperature)
            UB = [1, 1, Inf];
        elseif isequal(expMode, 'key')
            LB = [0, 0, -Inf]; % gamma (discount factor), alpha (learning rate), L (decision criterion)
            UB = [1, 1, Inf];
        end
        
    elseif isequal(Midx, 'SR_AttMod')
        theModel = @SR_model_choice_AttMod;
        if isequal(expMode, 'mouse')
            LB = [0, 0, 0, 0]; % gamma (discount factor), alpha (learning rate), eta (weight to update the learning rate trial-by-trial), beta (inverse temperature)
            UB = [1, 1, 1, Inf];
        end

    elseif isequal(Midx, 'SR_AttMod_adaptiveLR')
        theModel = @SR_model_choice_AttMod_adaptiveLR;
        if isequal(expMode, 'mouse')
            LB = [0, 0, 0, 0, 0]; % gamma (discount factor), alpha (learning rate), threshold, scaleFactor, beta (inverse temperature)
            UB = [1, 1, Inf, 1, Inf];
        end

    elseif isequal(Midx, 'SR_AttMod_adaptiveLR_nonLin')
        theModel = @SR_model_choice_AttMod_adaptiveLR_nonLin;
        if isequal(expMode, 'mouse')
            LB = [0, 0, 0, 0]; % gamma (discount factor), alpha (learning rate), eta (the exponential parameter), beta (inverse temperature)
            UB = [1, 1, Inf, Inf];
        end

    elseif isequal(Midx, 'SR_AttMod_lapse')
        theModel = @SR_model_choice_AttMod_lapse;
        if isequal(expMode, 'mouse')
            LB = [0, 0, 0, 0, 0]; % gamma (discount factor), alpha (learning rate), eta (weight to update the learning rate trial-by-trial), lambda (weight of random choice), beta (inverse temperature)
            UB = [1, 1, 1, 1, Inf];
        end

    elseif isequal(Midx, 'SRnet') % 'SRnet' are exactly the same as 'SR', the only difference is 'SRnet' only model a sequence with two options in each trial
        theModel = @SR_model_choice_networkLearnability;
        LB = [0, 0, 0]; % gamma (discount factor), alpha (learning rate), beta (inverse temperature)
        UB = [1, 1, Inf];
        
    elseif isequal(Midx, 'MaxEntropyNet') % 'MaxEntropyNet' are exactly the same as 'MaxEntropy'
        theModel = @maxEntropy_model_choice_networkLearnability;
        LB = [0, 0]; % beta (inverse temperature), r1 (inverse temperateure of the choice model)
        UB = [Inf, Inf];
        if isequal(priorWord, 'InformPrior')
            LB = [LB, 0]; % the initial value of the matrix nij
            UB = [UB, Inf];
        end

    elseif isequal(Midx, 'SRgd')  %% dynamic weight for the item matrix
        theModel = @SR_model_choice_gd;
        if isequal(expMode, 'mouse')
            LB = [0, 0, 0, 0, 0, 0]; 
            UB = [1, 1, 1, 1, 1, Inf];
        elseif isequal(expMode, 'key')
            LB = [0, 0, 0, 0, 0, -Inf]; 
            UB = [1, 1, 1, 1, 1, Inf];
        end
        
    elseif isequal(Midx, 'SRcw') %% constant weight for the item matrix
        theModel = @SR_model_choice_fixW;
        LB = [0, 0, 0, 0, 0, 0];
        UB = [1, 1, 1, 1, 1, Inf];
   
    elseif isequal(Midx, 'SRgdC')   
        theModel = @SR_model_choice_gdChoice;
        if isequal(expMode, 'mouse')
            LB = [0, 0, 0, 0, 0, 0, 0]; 
            UB = [1, 1, 1, 1, 1, Inf, Inf];
        elseif isequal(expMode, 'key')
            LB = [0, 0, 0, 0, 0, -Inf, -Inf]; 
            UB = [1, 1, 1, 1, 1, Inf, Inf];
        end
        
    elseif isequal(Midx, 'SRgdClusterTrace')   
        theModel = @clusterTrace_item_DynamicCompModel;
        LB = [0, 0, 0, 0, 0, 0];
        UB = [1, 1, Inf, Inf, 1, Inf];  
        
    elseif isequal(Midx, 'SRgdClusterWin') 
        theModel = @clusterWin_item_DynamicCompModel;
        LB = [0, 0, 2, 0, 0];
        UB = [1, 1, 200, 1, Inf];  
        
    elseif isequal(Midx, 'SRgdClusterWinDisc') 
        theModel = @clusterWinDisc_item_DynamicCompModel;
        LB = [0, 0, 2, 0, 0, 0];
        UB = [1, 1, 200, 1, 1, Inf];

    elseif isequal(Midx, 'SRsymm')
        theModel = @SRsymm_model_choice;
        if isequal(expMode, 'key')
            LB = [0, 0, -Inf]; % gamma (discount factor), alpha (learning rate), L (decision criterion)
            UB = [1, 1, Inf];
        end
        
    elseif isequal(Midx, 'heurProcess')
        theModel = @heurProcess_model_choice;
        if isequal(expMode, 'mouse')
            LB = [2, 0, 0]; % winLen (time window length), inverse temperature, alpha_mw
            UB = [200, Inf, 1];
        elseif isequal(expMode, 'key')
            LB = [2, -Inf]; % winLen (time window length), L (decision criterion)
            UB = [200, Inf];
        end
        
    elseif isequal(Midx, 'heurProcessProb')
        theModel = @heurProcessProb_model_choice;
        if isequal(expMode, 'mouse')
            LB = [2, 0, 0]; % winLen (time window length), inverse temperature, alpha_mw
            UB = [200, Inf, 1];
        end
        
    elseif isequal(Midx, 'heurProcessDiscount')
        theModel = @heurProcessDiscount_model_choice;
        if isequal(expMode, 'mouse')
            LB = [2, 0, 0, 0]; % winLen (time window length), gamma (discounting), inverse temperature, alpha_mw
            UB = [200, 1, Inf, 1];
        elseif isequal(expMode, 'key')
            LB = [2, 0, -Inf]; % winLen (time window length), gamma (discounting), L (decision criterion)
            UB = [200, 1, Inf];
        end
        
    elseif isequal(Midx, 'heurProcessTrace')
        theModel = @heurProcessTrace_model_choice;
        if isequal(expMode, 'mouse')
            LB = [0, 0, 0];  % beta_cls, beta, alpha_mw
            UB = [Inf, Inf, 1];
        end

    elseif isequal(Midx, 'heurProcessTraceDecay')
        theModel = @heurProcessTrace_model_choice_withDecay;
        if isequal(expMode, 'mouse')
            LB = [0, 0, 0];  % beta_cls, beta, alpha_mw
            UB = [Inf, Inf, 1];
        end

    elseif isequal(Midx, 'heurProcessTrace_AttMod')
        theModel = @heurProcessTrace_model_choice_AttMod;
        if isequal(expMode, 'mouse')
            LB = [0, 0, 0, 0];  % beta_cls, beta, eta (weight to update the learning rate trial-by-trial), alpha_mw
            UB = [Inf, Inf, 1, 1];
        end

    elseif isequal(Midx, 'heurProcessTrace_AttMod_adaptiveLR')
        theModel = @heurProcessTrace_model_choice_adaptiveLR;
        if isequal(expMode, 'mouse')
            LB = [0, 0, 0, 0, 0];  % beta_cls, beta, alpha_mw, threshold, scaleFactor
            UB = [Inf, Inf, 1, Inf, 1];
        end

    elseif isequal(Midx, 'heurProcessTrace_AttMod_adaptiveLR_nonLin')
        theModel = @heurProcessTrace_model_choice_adaptiveLR_nonLin;
        if isequal(expMode, 'mouse')
            LB = [0, 0, 0, 0];  % beta_cls, beta, alpha_mw, eta (the exponential parameter)
            UB = [Inf, Inf, 1, Inf];
        end

    elseif isequal(Midx, 'heurProcessTraceAvg')
        theModel = @heurProcessTraceAvg_model_choice;
        if isequal(expMode, 'mouse')
            LB = [0, 0, 0];  % beta_cls, beta, initV (the initial value of the matrix nij)
            UB = [Inf, Inf, Inf];
        end
        
    elseif isequal(Midx, 'heurProcessTraceWin')
        theModel = @heurProcessTraceWin_model_choice;
        if isequal(expMode, 'mouse')
            LB = [2, 0, 0];  % winLen (time window length), beta, alpha_mw
            UB = [200, Inf, 1];
        end
  
    elseif isequal(Midx, 'SRheurProcessProb') 
        theModel = @SR_heurProcessProb_choice_gdChoice;
        if isequal(expMode, 'mouse')
            LB = [0, 0, 2, 0, 0, 0, 0]; % gamma_trans (discount factor), alpha_trans (learning rate), winLen (time window length), alpha_mw (weight of the delta-rule updating), alpha_gd (learning rate of gradient descent), beta_trans, beta_ass
            UB = [1, 1, 200, 1, 1, Inf, Inf];
        end
        
    elseif isequal(Midx, 'SRheurProcessTrace') %%% ******** The best model ********
        theModel = @SR_heurProcessTrace_choice_gdChoice;
        if isequal(expMode, 'mouse')
            LB = [0, 0, 0, 0, 0, 0, 0]; % gamma_trans (discount factor), alpha_trans (learning rate), beta_cls (inverse temperature), alpha_mw (weight of the delta-rule updating), alpha_gd (learning rate of gradient descent), beta_trans, beta_ass
            UB = [1, 1, Inf, 1, 1, Inf, Inf];
        end

    elseif isequal(Midx, 'SRheurProcessTraceTrunc')
        theModel = @SR_heurProcessTrace_choice_gdChoice_trunc;
        if isequal(expMode, 'mouse')
            LB = [0, 0, 0, 0, 0, 0, 0]; % gamma_trans (discount factor), alpha_trans (learning rate), beta_cls (inverse temperature), alpha_mw (weight of the delta-rule updating), alpha_gd (learning rate of gradient descent), beta_trans, beta_ass
            UB = [1, 1, Inf, 1, 1, Inf, Inf];
        end
        
    elseif isequal(Midx, 'SRheurProcessTraceTruncAttmod')
        theModel = @SR_heurProcessTrace_choice_gdChoice_truncAttMod;
        if isequal(expMode, 'mouse')
            LB = [0, 0, 0, 0, 0, 0, 0, 0]; % gamma_trans (discount factor), alpha_trans (learning rate), beta_cls (inverse temperature), alpha_mw (weight of the delta-rule updating), alpha_gd (learning rate of gradient descent), eta (weight to update the learning rate trial-by-trial), beta_trans, beta_ass
            UB = [1, 1, Inf, 1, 1, 1, Inf, Inf];
        end

    elseif isequal(Midx, 'SRheurProcessTraceTruncAtt')
        theModel = @SR_heurProcessTrace_choice_gdChoice_truncAtt;
        if isequal(expMode, 'mouse')
            LB = [0, 0, 0, 0, 0, 0, 0]; % gamma_trans (discount factor), alpha_trans (learning rate), beta_cls (inverse temperature), alpha_mw (weight of the delta-rule updating), alpha_gd (learning rate of gradient descent), beta_trans, beta_ass
            UB = [1, 1, Inf, 1, 1, Inf, Inf];
        end

    elseif isequal(Midx, 'SRheurProcessTraceTruncAttDyn') % dynamically modulate the updating parameter in the association learning component
        theModel = @SR_heurProcessTrace_choice_gdChoice_truncAttDyn;
        if isequal(expMode, 'mouse')
            LB = [0, 0, 0, 0, 0, 0, 0]; 
            UB = [1, 1, Inf, 1, 1, Inf, Inf];
        end

    elseif isequal(Midx, 'SRheurProcessTraceTruncAttFilter') 
        theModel = @SR_heurProcessTrace_choice_gdChoice_truncAttFilter;
        if isequal(expMode, 'mouse')
            LB = [0, 0, 0, 0, 0, 0, 0, -Inf];
            UB = [1, 1, Inf, 1, 1, Inf, Inf, Inf];
        end

    elseif isequal(Midx, 'SRheurProcessTraceExpDecay') % same parameters as 'SRheurProcessTraceTruncAttDyn'
        theModel = @SR_heurProcessTrace_choice_gdChoice_expDecay;
        if isequal(expMode, 'mouse')
            LB = [0, 0, 0, 0, 0, 0, 0, 0]; % last parameter: d_ass (the decay parameter)
            UB = [1, 1, Inf, 1, 1, Inf, Inf, Inf];
        end

    elseif isequal(Midx, 'SRheurProcessTraceLinDecay') % linear decay
        theModel = @SR_heurProcessTrace_choice_gdChoice_linDecay;
        LB = [0, 0, 0, 0, 0, 0]; % gamma_trans (discount factor), alpha_trans (learning rate), k (linear decay parameter), alpha_gd (learning rate of gradient descent), beta_trans, beta_ass
        UB = [1, 1, 1, 1, Inf, Inf];

    elseif isequal(Midx, 'SRheurProcessTraceCluster')
        theModel = @SR_heurProcessTrace_choice_gdChoice_cls; % its parameter space is exactly the same as SR_heurProcessTrace_choice_gdChoice model
        if isequal(expMode, 'mouse')
            LB = [0, 0, 0, 0, 0, 0, 0]; 
            UB = [1, 1, Inf, 1, 1, Inf, Inf];
        end

    elseif isequal(Midx, 'SRheurProcessTraceAvg')
        theModel = @SR_heurProcessTraceAvg_choice_gdChoice;
        if isequal(expMode, 'mouse')
            LB = [0, 0, 0, 0, 0, 0, 0]; % gamma_trans (discount factor), alpha_trans (learning rate), beta_cls (inverse temperature), initV (the initial value of the matrix nij), alpha_gd (learning rate of gradient descent), beta_trans, beta_ass
            UB = [1, 1, Inf, Inf, 1, Inf, Inf];
        end

    elseif  isequal(Midx, 'SRheurProcessTraceWin')
        theModel = @SR_heurProcessTraceWin_choice_gdChoice;
        if isequal(expMode, 'mouse')
            LB = [0, 0, 2, 0, 0, 0, 0]; % gamma_trans (discount factor), alpha_trans (learning rate), winLen (time window length), alpha_mw (weight of the delta-rule updating), alpha_gd (learning rate of gradient descent), beta_trans, beta_ass
            UB = [1, 1, 200, 1, 1, Inf, Inf];
        end

    elseif isequal(Midx, 'SRheurProcessTraceDistort')
        theModel = @SR_heurProcessTraceDistort_choice_gdChoice;
        if isequal(expMode, 'mouse')
            LB = [0, 0, 0, 0, 0, 0, 0, 0]; % gamma_trans (discount factor), alpha_trans (learning rate), beta_cls (inverse temperature), alpha_mw (weight of the delta-rule updating), alpha_gd (learning rate of gradient descent), beta_trans, beta_ass, r (exponential term of the one-parameter Prelec)
            UB = [1, 1, Inf, 1, 1, Inf, Inf, 1];
        end
        
    elseif isequal(Midx, 'objTP')
        theModel = @objTP_model_choice;
        if isequal(expMode, 'mouse')
            LB = [0]; % inverse temperature
            UB = [Inf];
        elseif isequal(expMode, 'key')
            LB = [-Inf]; % L (decision criterion)
            UB = [Inf];
        end
        
    elseif isequal(Midx, 'heurTP')
        theModel = @heurTP_model_choice;
        if isequal(expMode, 'key')
            LB = [1e-16, 1e-16, -Inf]; % L (decision criterion)
            UB = [5, 5, Inf];
        end
    %%%========== 3 component model ==========
    %SR_Ass_Oneshot_choice_gdChoice_trunc
    elseif isequal(Midx, 'SRAssOneshotLin3Comp')
        theModel = @SR_Ass_Oneshot_choice_gdChoice_trunc;
        if isequal(expMode, 'mouse')
            LB = [0, 0, 0, 0, 0, 0, 0, 0, 0];
            UB = [1, 1, Inf, 1, 1, 1, Inf, Inf, Inf];
        end
         % gamma_trans (discount factor), alpha_trans (learning rate), 
         % beta_cls (inverse temperature), alpha_mw (weight of the delta-rule updating), 
         % k (linear decay parameter in the one-shot learning component)
         % alpha_gd (learning rate of gradient descent)
         % eta: parameter to control the weight of one-shot learning
         % beta_trans, beta_ass

    %%%========== one-shot learning related ==========
    elseif isequal(Midx, 'oneShotLinDecay')
        theModel = @oneShot_linearDecay;
        if isequal(expMode, 'mouse')
            LB = [0, 0]; % beta (inverse temperature), k (linear decay parameter)
            UB = [Inf, 1];
        end

    elseif isequal(Midx, 'oneShotExpDecay')
        theModel = @oneShot_expDecay;
        if isequal(expMode, 'mouse')
            LB = [0, 0]; % beta (inverse temperature), d_ass (the decay parameter)
            UB = [Inf, Inf];
        end

    elseif isequal(Midx, 'SRoneShotLinDecay') 
        theModel = @SR_oneShot_linearDecay;
        if isequal(expMode, 'mouse')
            LB = [0, 0, 0, 0, 0, 0]; % gamma_trans (discount factor), alpha_trans (learning rate), k (linear decay parameter in the one-shot learning component), alpha_gd (learning rate of gradient descent), beta_trans, beta_ass
            UB = [1, 1, 1, 1, Inf, Inf];
        end

    elseif isequal(Midx, 'SRoneShotExpDecay')
        theModel = @SR_oneShot_expDecay;
        if isequal(expMode, 'mouse')
            LB = [0, 0, 0, 0, 0, 0]; % gamma_trans (discount factor), alpha_trans (learning rate), d_ass (decay parameter of the one-shot learning), alpha_gd (learning rate of gradient descent), beta_trans, beta_ass
            UB = [1, 1, Inf, 1, Inf, Inf];
        end

    end
    options = optimset('MaxFunEvals', length(LB)*1000, 'MaxIter', length(LB)*1000);
    %options = optimset('MaxFunEvals', length(LB)*10000, 'MaxIter', length(LB)*10000);
%     if isequal(Midx, 'SRgd') || isequal(Midx, 'SRgdC')
%         options = optimset('MaxFunEvals', length(LB)*1000, 'MaxIter', length(LB)*1000);
%     else
%         options = optimset('MaxFunEvals', length(LB)*40000, 'MaxIter', length(LB)*40000);
%     end
    iFit    = 1;
    lli_all    = zeros(1, nFit);
    params_all = zeros(length(LB), nFit);
    while iFit <= nFit
        display([expID, '-', subID, '-', Midx, 'model-', fitWord, '-Fit', num2str(iFit)]);
        if isequal(Midx, 'heurProcess') || isequal(Midx, 'heurProcessProb') || isequal(Midx, 'heurProcessDiscount') || isequal(Midx, 'heurProcessTraceWin')
            x0 = [ceil(1+10*rand), rand(1, length(LB)-1)];
        else
            x0 = rand(1, length(LB));
        end
        tic
        [paramsEst, minuslli, exitflag] = fminsearchbnd(@(params)theModel(params, stim, resp, states, choiceId, delta_t_trials, assDist_mat, expMode, posXY, respYes, mouseTraj, objAng_col, objDtrNo_col, miniDmat, fitWord, transMat_input, nij_input, priorWord), x0, LB, UB, options)
        toc
        
        if exitflag == 1
            params_all(:, iFit) = paramsEst;
            lli_all(iFit)       = - minuslli;
            iFit = iFit + 1;
        end
    end
    [bestlli, llibest_id] = max(lli_all);
    
    paramsEst = params_all(:, llibest_id);
    lliEst    = bestlli;
    save(fn, 'lliEst', 'paramsEst', 'lli_all');
end


