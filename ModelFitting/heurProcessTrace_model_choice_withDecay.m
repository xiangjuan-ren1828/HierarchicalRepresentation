function [minuslli, M_str, p_choice_all, choiceID_all] = heurProcessTrace_model_choice_withDecay(params, stim, resp, states, choiceId, delta_t_trials, assDist_mat, expMode, posXY, respYes, mouseTraj, objAng, objDtrNo, miniDmat, fitWord, transMat_input, nij_input, priorWord)
% write by XR @ Jan 13 2026
% modified based on heurProcessTrace_model_choice.m: for the non-displayed
% association pairs on the current trial, decay will be implemented
% =========================================================================
% write by rxj @ 06/29/2022 based on heurProcess_model_choice.m
% Using the same idea as Lynn et al., 2020, but assuming a symmetrical
% transition matrix
% =========================================================================
% write by rxj @ 08/09/2021 based on the followings
% function [] = SR_model()
% =========================================================================
% write by rxj @ 05/05/2021 based on the followings
% function SR_fit().m
% write by rxj @ 07/24/2020
% SR fitting
%
% stim: from nodes
% resp: to nodes; for the keypress experiment, it must be 2 columns

if isequal(expMode, 'mouse')
    %%% parameters for item association update
    beta_cls  = params(1); % inverse temperature of maxEntropy model
    %%% parameters for probabilistic choice
    beta      = params(2); % inverse temperature
    Z         = 1 / (1 - exp(-beta_cls));
    %initV = params(3);
    %nij_input(1 : (states + 1) : end) = initV;
    alpha_mw = params(3); % weight of the delta-rule updating
    
    %%
    trialLen = length(stim);
    %nij      = nij_input;      %zeros(states, states);
    M_str    = transMat_input; %eye(states, states);
    minuslli      = 0;
    choiceID_all = nan(trialLen, max(objDtrNo) + 1);
    p_choice_all = nan(trialLen, max(objDtrNo) + 1);
    %%
    for iT = 1 : trialLen
        s_i    = stim(iT);    % from which node
        s_Tgt  = resp(iT, 1); % to which node
        s_nTgt = resp(iT, 2);
        
        %%% index of variable
        fromNd_i    = stim(iT);
        toNd_i      = resp(iT, :);
        posX        = posXY(iT, :, 1);
        posY        = posXY(iT, :, 2);
        ang_i       = objAng(iT, :);       %% target + distractor
        dtrNo_i     = objDtrNo(iT, :) + 1; %% number of object (including target)
        %%% index according to distractor numbers
        toNd_i      = toNd_i(1 : dtrNo_i);
        posX        = posX(1 : dtrNo_i);
        posY        = posY(1 : dtrNo_i);
        ang_i       = ang_i(1 : dtrNo_i);
        respYes_i   = respYes{iT};
        mouseTraj_i = mouseTraj{iT};
        choiceRorW  = choiceId(iT);
        [rankID, dtrChoiceID] = multiObj_choice_match(fromNd_i, toNd_i, posX, posY, respYes_i, mouseTraj_i, ang_i, dtrNo_i, choiceRorW, miniDmat, fitWord);
        
        %% item association learning: directional
        %%% assuming the matrix is symmetrical
        Vopt = nan(dtrNo_i, 1);
        for jdtr = 1 : dtrNo_i
            id_Tgt = sort([s_i, toNd_i(jdtr)]);
            Vopt(jdtr) = M_str(id_Tgt(2), id_Tgt(1));
        end
        p_denominator = exp(beta * Vopt);
        p_denominator(isinf(p_denominator)) = exp(700);
        pCorr = p_denominator ./ sum(p_denominator);
        pCorr(pCorr < 1e-16)   = 1e-16;
        pCorr(pCorr > 1-1e-16) = 1-1e-16;
        
        choiceID_all(iT, 1 : length(pCorr)) = dtrChoiceID;
        p_choice_all(iT, 1 : length(pCorr)) = pCorr;
        
        %% Clustering learning: nondirectional
        % update the 'Clustering'transition matrix
        Bt = zeros(states, 1);        
        halfMask    = tril(true(states), -1); % the half you store: strictly lower triangle
        updatedMask = false(states, states);  % mark which pairs got updated this trial
        for iS = 1 : states
            delta_t = delta_t_trials{iT, iS};
            if ~isempty(delta_t)
                P_delta = exp(-beta_cls .* delta_t) ./ Z;
                Bt(iS)  = sum(P_delta);
                id_Tgt  = sort([iS, s_Tgt]);
                % ------ update the current association pair ------
                M_str(id_Tgt(2), id_Tgt(1)) = (1 - alpha_mw) * M_str(id_Tgt(2), id_Tgt(1)) + alpha_mw * Bt(iS);

                % ------ mark as updated ------
                updatedMask(id_Tgt(2), id_Tgt(1)) = true;
            end
        end  
        % ---- decay once for all other pairs in the same half ----
        decayMask = halfMask & ~updatedMask;
        M_str(decayMask) = (1 - alpha_mw) * M_str(decayMask);
        
        %% calculate the likelihood
        minuslli = minuslli + sum(dtrChoiceID .* log(pCorr));
        
    end
    minuslli = -minuslli;

elseif isequal(expMode, 'key')
    winLen = params(1);  % time window length, integer & no less than 2
    winLen = round(winLen); 
    %gamma  = params(2);
    L      = params(2);  % the decision criterion 
    sigma  = 1;
    
    %%
    trialLen = length(stim);
    M_str  = transMat_input; 
    heurTransCell = cell(states, states);
    minuslli      = 0;
    choiceID_all  = choiceId;
    p_choice_all  = nan(trialLen, 3);
    %%
    for iT = 1 : trialLen
        s_i    = stim(iT);
        s_Tgt  = resp(iT, 1);
        s_nTgt = resp(iT, 2);
        %V_Tgt  = M(s_i, s_Tgt);
        %V_nTgt = M(s_i, s_nTgt);
        
        %%% assuming the matrix is symmetrical
        id_Tgt  = sort([s_i, s_Tgt]);
        id_nTgt = sort([s_i, s_nTgt]);
        V_Tgt   = M_str(id_Tgt(2), id_Tgt(1));
        V_nTgt  = M_str(id_nTgt(2), id_nTgt(1));

        %% ===============================================================
        %%% 3 kinds of choice (alike signal detection theory): (1) if
        %%% |Vt-Vd|<=L, no response; (2) if Vt-Vd<-L, incorrect response
        %%% before cue; (3) if Vt-Vd>L, correct response before cue
        %%% ===============================================================
        p_choice = zeros(3, 1);
        Vd = V_Tgt-V_nTgt;
        mu = Vd;
        pd = makedist('Normal', 'mu', mu, 'sigma', sigma);
        % plot: x = -5:0.1:5; y = pdf(pd, x); plot(x, y);
        %L1 = Lraw1 * abs(Vd);
        %L2 = Lraw2 * abs(Vd);

        p_choice(2) = cdf(pd, 0-L);     % ..........incorrect response before cue (selecting the non-target)
        p_choice(3) = 1 - cdf(pd, 0+L); % ..........correct response before cue (selecting the target)
        p_choice(1) = cdf(pd, 0+L) - cdf(pd, 0-L); % response after cue
        dtrChoiceID = choiceId(iT, :); 
        %%% ===the 'choiceId' in expMode='mouse' & expMode='key' are different===
        p_choice_all(iT, :) = p_choice;
        
        %% update the successor representation matrix
        if iT < winLen
            tgtWin_i = resp(1 : iT, 1);
        else
            tgtWin_i = resp((iT - winLen + 1) : iT, 1);
        end
        
        %% only with a parameter: time window length
        s_uniq = unique(tgtWin_i);
        if length(s_uniq) == 1
            s_uniq = [stim(iT); s_uniq];
        end
        s_pair = nchoosek(s_uniq, 2);
        for ii = 1 : size(s_pair, 1)
            iFr = s_pair(ii, 1);
            iTo = s_pair(ii, 2);
            %%% item frequency
            itemFreq = length(find(tgtWin_i == iFr | tgtWin_i == iTo)); %% frequency of individual items
            
            %%% pair frequency
            idx_mark = nan(length(tgtWin_i), 1);
            idx_mark((tgtWin_i == iFr)) = 1;
            idx_mark((tgtWin_i == iTo)) = 2;
            idx_mark(isnan(idx_mark)) = [];
            
            pairFreq = length(find((diff(idx_mark)) ~= 0));
            heurTransCell{iTo, iFr} = [heurTransCell{iTo, iFr}; pairFreq/itemFreq];
            
            %%% running sum across trials
            pii = heurTransCell{iTo, iFr};
            if isempty(pii)
                M_str(iTo, iFr) = 0;
            else
                M_str(iTo, iFr) = nansum(pii); %nanmean(heurTransCell{iTo, iFr, iL});
            end
        end
        %%% normalization
        if (max(M_str(:)) - min(M_str(:))) ~= 0
            M_str = (M_str - min(M_str(:))) ./ (max(M_str(:)) - min(M_str(:))); % rescaling data between 0-1
        end
        
        %% calculate the likelihood
        if p_choice(dtrChoiceID == 1) < 1e-16
            p_choice(dtrChoiceID == 1) = 1e-16;
        elseif p_choice(dtrChoiceID == 1) > 1-1e-16
            p_choice(dtrChoiceID == 1) = 1-1e-16;
        end
        minuslli = minuslli + dtrChoiceID(dtrChoiceID == 1) .* log(p_choice(dtrChoiceID == 1));
        
    end
    minuslli = -minuslli;
    % the Ncorr from the trud data was generated according to a bimomial distribution, that is, Ncorr=binornd(Nall, pCorr)
    %minuslli = -sum(choiceId .* log(pCorr) + (1 - choiceId) .* log(1 - pCorr)); % the likelihood function describe the relatinship between true data and your model
end
