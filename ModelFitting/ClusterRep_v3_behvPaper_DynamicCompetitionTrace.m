% ClusterRep_v3_behvPaper_DynamicCompetitionTrace.m
% modified by XR @ Jan 27 based on the followings
% ========================================================================
% DynamicCompetitionTrace.m
% write by rxj @ 06/07/2023
% For the mixture model of cluster learning and item learning, plot the
% dynamic change of the weight

clear
clc

%%
addpath('tight_subplot/');
addpath(genpath('HierarchicalCluster/'));

%% Columns index for the clusterResult.mat
rndTrial = 700;
HamTrial = 800;
nTrials  = rndTrial + HamTrial;
nBin     = 100;
BinL     = nTrials / nBin;

expModeList = {'mouse', 'key'};
iMode       = 1;
expMode_i   = expModeList{iMode};
if isequal(expMode_i, 'mouse')
    rndOrHam = 1; %% trial countings for the random (1-700) and hamiltonian (1-800) walk
    blockNo  = 2; %% block number
    trialNo  = 3; %% trial number
    trialTp  = 4; %% trial type: 1-random; 2-hamiltonian
    objTgt   = 5; %% target object
    objResp  = 6 : 9; %% responded objects
    moveTime = 10; %% moving time of the mouse
    respRT   = 11; %% response reaction time
    dispRT   = 12; %% object total display time
    errorIdx = 13; %% whether error (including incorrect response and time out) or not: 1-yes
    timeOIdx = 14; %% whether time out: 1-yes
    bonusCn  = 15; %% bonus in each trial
    trialLt  = 16; %% trial total length
    imgPosX  = 17 : 20;
    imgPosY  = 21 : 24;
    fastIdx  = 25;
elseif isequal(expMode_i, 'key')
    rndOrHam = 1; %% trial countings for the random (1-700) and hamiltonian (1-800) walk
    blockNo  = 2; %% block number
    trialNo  = 3; %% trial number
    trialTp  = 4; %% trial type: 1-random; 2-hamiltonian
    objTgt   = 5; %% target object
    objResp  = 6 : 7; %% responded objects
    cueOn    = 8; %% cue onset time after stimuli onsets: jitter between 0.9 and 1.2
    duraT    = 9;
    moveTime = 10; %% moving time of the mouse
    respRT   = 11; %% response reaction time
    dispRT   = 12; %% object total display time
    errorIdx = 13; %% whether error (including incorrect response and time out) or not: 1-yes
    timeOIdx = 14; %% whether time out: 1-yes
    errCnt   = 15;
    bonusCn  = 16; %% bonus in each trial
    trialLt  = 17; %% trial total length
    imgPosX  = 18 : 19;
    imgPosY  = 20 : 21;
    fastIdx  = 22;
    imgAngs  = 23 : 24;
    fixLast  = 25;
    base_err = 26 : 27;
    resp_err = 28 : 29;
    fixFroze = 30;
    keyPress = 31;
end
miniDmat = modelRDM_minimalDist();
clsDef  = [1, 2, 3, 4, 5; ...  %% the 1st cluster
           6, 7, 8, 9, 10; ... %% the 2nd cluster
           11, 12, 13, 14, 15];%% the 3rd cluster

%% Extracting relevant parameters from the model
% added by rxj @ 08/12/2021
fitWord   = 'angleChoice'; % {'clickChoice', 'angleChoice'}
iteWord   = 'naive';       % {'naive', 'trained'}
trialWord = 'all';         % {'all', 'Hamiltonian'}
states = 15;
nFit   = 100;
expList   = {'ImplicitExp', 'ExplicitExp', 'ImplicitRandExp'};
MList     = {'SRheurProcessTraceTruncAttMod'}; % 'SRheurProcessProb', 'SRheurProcessTraceWin', 'SRheurProcessTrace', 'SRheurProcessTraceTrunc'
priorList = {'nullPrior'}; % 'nullPrior', 'nullPrior', 'nullPrior'
folder    = '/Users/ren/Projects-NeuroCode/MyExperiment/HierarchicalCluster';
refit     = 0;
calFlg    = 0; % 0: one-time; 1: accumulated
nanRemove = 1; % 1: remove 'nan' from choiceId and other variales before importing them into fitting procedure
nSim = 500;

%%
subLen = 24;
lambda_trials_subj_nExp     = nan(subLen, nTrials, length(expList));
lambda_mean_subj_nExp       = nan(subLen, 4, length(expList)); % Col 1-2: within and between transition across all trials; Col3-4: only the last 800 trials with Random and Hamiltonian walks intermixed
lambda_mean_trans_walk_nExp = nan(subLen, 2, 2, length(expList)); % (only the last 800 trials) first 2: within vs. between trans; second 2: random vs. ham walk
for iExp = 1 : length(expList)
    ExpWord = expList{iExp};
    %% Subject
    if isequal(ExpWord, 'ImplicitExp')
        bhvDataDir = [folder, '/FormalExp-LynnNetwork-Results/'];
        subj_list  = {'wsn_1_f_18', 'dy_2_f_22', 'haq_3_f_24', 'hry_4_f_20', 'zjx_5_m_20', 'yyq_6_f_18', 'zkw_7_m_18', 'zy_8_f_20', 'hys_9_m_20', 'cjj_10_m_18', ...
                      'dwq_11_f_22', 'ljl_12_m_20', 'jyx_13_m_19', 'zk_14_f_21', 'lsy_15_m_19', 'cjl_16_m_19', 'yjy_17_f_23', 'lym_18_f_19', 'pr_19_f_23', 'ws_20_f_21', ...
                      'wn_21_f_21', 'hjy_22_f_18', 'qyk_23_f_22', 'yd_24_f_20'};
        subjLab    = {'wsn1', 'dy2', 'haq3', 'hry4', 'zjx5', 'yyq6', 'zkw7', 'zy8', 'hys9', 'cjj10', ...
                      'dwq11', 'ljl12', 'jyx13', 'zk14', 'lsy15', 'cjl16', 'yjy17', 'lym18', 'pr19', 'ws20', ...
                      'wn21', 'hjy22', 'qyk23', 'yd24'};
        nBlock = 7;
        
    elseif isequal(ExpWord, 'ExplicitExp')
        bhvDataDir = [folder, '/FormalExp-ExplicitLearning-LynnNetwork-Results/'];
        subj_list  = {'hsp_1_m_21', 'hr_2_f_22', 'pxy_3_m_19', 'wyx_4_m_22', 'lml_5_f_19', 'lf_6_m_21', 'md_7_f_21', 'srz_8_m_20', 'fsq_9_f_18', 'caq_10_f_20', ...
                      'xjm_11_f_26', 'cjs_12_m_19', 'xxx_13_f_19', 'wsq_14_f_19', 'zzm_15_f_18', 'lsy_16_f_20', 'man_17_f_19', 'zxy_18_f_18', 'cyh_19_f_21', 'szn_20_f_19', ...
                      'sb_21_m_22', 'drq_22_f_18', 'cy_23_f_18', 'zhc_24_m_21'};
        subjLab    = {'hsp1', 'hr2', 'pxy3', 'wyx4', 'lml5', 'lf6', 'md7', 'srz8', 'fsq9', 'caq10', ...
                      'xjm11', 'cjs12', 'xxx13', 'wsq14', 'zzm15', 'lsy16', 'man17', 'zxy18', 'cyh19', 'szn20', ...
                      'sb21', 'drq22', 'cy23', 'zhc24'};
        nBlock = 7; 
        
    elseif isequal(ExpWord, 'ImplicitRandExp')
        bhvDataDir = [folder, '/FormalExp-LynnNetwork-ImplicitRandom-Results/'];
        subj_list  = {'zyh_1_f_21', 'why_2_m_18', 'cr_3_f_19', 'zyx_4_f_25', 'wym_5_m_23', 'wd_6_f_18', 'lyh_7_f_21', 'zr_8_f_24', 'zyh_9_m_19', 'zzy_10_f_19', ...
                      'smq_11_f_25', 'sz_12_f_20', 'lzy_13_m_22', 'yxy_14_f_18', 'zxl_15_f_18', 'wqh_16_m_23', 'zxj_17_m_21', 'skx_18_f_20', 'zlh_19_m_24', 'gwt_20_f_23', ...
                      'lwn_21_f_18', 'lrp_22_m_21', 'sjj_23_f_18', 'xy_24_f_19'}; % , 'xr_18_f_23'
        subjLab    = {'zyh1', 'why2', 'cr3', 'zyx4', 'wym5', 'wd6', 'lyh7', 'zr8', 'lyh9', 'zzy10', ...
                      'smq11', 'sz12', 'lzy13', 'yxy14', 'zxl15', 'wqh16', 'zxj17', 'skx18', 'zlh19', 'gwt20', ...
                      'lwn21', 'lrp22', 'sjj23', 'xy24'}; % , 'xr18'
        nBlock = 7;

    end
    if isequal(expMode_i, 'mouse') %% read the mouseTraj.csv
        angleDir = [bhvDataDir, 'angleCal/'];
        angleTraj_pd_subj = readtable([angleDir, 'angleTraj_pd_subj.csv']);
        %%% variables
        subNo_col     = angleTraj_pd_subj.subNo;
        trialCnt_col  = angleTraj_pd_subj.trialCnt;
        blockNo_col   = angleTraj_pd_subj.blockNo;
        trialNo_col   = angleTraj_pd_subj.trialNo;
        tgtAng_col    = angleTraj_pd_subj.tgtAng;
        dtrAng1_col   = angleTraj_pd_subj.dtrAng1;
        dtrAng2_col   = angleTraj_pd_subj.dtrAng2;
        dtrAng3_col   = angleTraj_pd_subj.dtrAng3;
        dtrNoCnt_col  = angleTraj_pd_subj.dtrNoCnt; %% 1, 2, 3
    end
    subj_listBv = subj_list;
    subLen      = length(subj_list);
    
    %% Loop over participants
    M_trans_subj       = nan(15, 15, subLen, length(MList));
    M_ass_subj         = nan(15, 15, subLen, length(MList));

    %hiddenMat = zeros(states, states, 2, subLen); % 2 means 2 models: maxEntropy and SR  
    hiddenMat       = cell(subLen, 2);
    transStyle_subj = cell(subLen, 1);
    rndOrHam_subj   = cell(subLen, 1); 
    dtNum_subj      = cell(subLen, 1);
    for iSub = 1 : subLen
        disp([ExpWord, '-subj', num2str(iSub)]);
        %%
        subID  = subjLab{iSub};
        subjBv = subj_listBv{iSub};
        subjDir = [bhvDataDir, subjBv, '/'];
        
        if isequal(expMode_i, 'mouse')
            load([subjDir, subjBv, 'clusterResult_Blocks.mat'], 'clusterResult');
            load([subjDir, subjBv,  'respYes_trials_blc.mat'], 'respYes_trials_blc');
            load([subjDir, subjBv,  'mouseTraj_trials_blc.mat'], 'mouseTraj_trials_blc');
            %%% !!!!!!!! remove the time-out trials !!!!!!!! %%%
            rndOrHam = clusterResult(:, trialTp); %% trial type: 1-random; 2-hamiltonian
            stim     = clusterResult(:, objTgt);
            resp     = clusterResult(:, objResp);
            posXY        = zeros(length(stim), length(imgPosX), 2);
            posXY(:, :, 1) = clusterResult(:, imgPosX);
            posXY(:, :, 2) = clusterResult(:, imgPosY);
            errorIdx_Col   = clusterResult(:, errorIdx); %% incorrect click or time-out
            timeOIdx_Col   = clusterResult(:, timeOIdx);
            fast_Col       = clusterResult(:, fastIdx);
            exId           = find(errorIdx_Col ~= 0);
            choiceId_click = ones(length(stim), 1);
            choiceId_click(exId) = 0;
            dtNum = arrayfun(@(x) length(find(resp(x, :)~=0)), (1 : length(stim))');

            %%% transitional type
            from_nodes = clusterResult(:, objTgt);
            to_nodes   = clusterResult(:, objResp(1));
            transIn    = zeros(nTrials, 3);
            for ic = 1 : size(clsDef, 1)
                if ic == 1
                    clsIdx = find(from_nodes >= 1 & from_nodes <= 5);
                elseif ic == 2
                    clsIdx = find(from_nodes >= 6 & from_nodes <= 10);
                elseif ic == 3
                    clsIdx = find(from_nodes >= 11 & from_nodes <= 15);
                end
                transIn(clsIdx, 1) = arrayfun(@(x) ismember(x, clsDef(ic, :)), from_nodes(clsIdx));
                transIn(clsIdx, 2) = arrayfun(@(x) ismember(x, clsDef(ic, :)), to_nodes(clsIdx));
                transIn(clsIdx, 3) = repmat(ic, length(clsIdx), 1);
            end
            transStyle = (transIn(:, 1) == transIn(:, 2));
            transStyle_history = nan(length(stim), 1); % the transitional style before current trial
            for jT = 2 : length(stim)
                transStyle_history(jT) = transStyle(jT - 1);
            end
   
            %%% prepare the delta_t for Maxmum Entropy Model, which will be
            %%% transmitted to the function maxEntropy_model_choice.m
            delta_t_trials = cell(length(stim), states);
            for iT = 1 : length(stim)
                s_hs = stim(1 : iT); %% stimulus history
                for iS = 1 : states
                    delta_t = find(s_hs == iS);
                    if ~isempty(delta_t)
                        delta_t = iT - delta_t;
                        delta_t_trials{iT, iS} = delta_t;
                    end
                end
            end

            %%% distance between association pair
            assDist_mat = assDist_rec_cal(stim, states, calFlg);

            %%% distance of the occurrence time of an association pair to
            %%% the current trial number
            [relation_mat, relation_history] = oneShot_relation_cal(stim, resp, states, miniDmat);
            
            choiceId_ang = nan(length(stim), 1);
            objAng_col   = zeros(length(stim), 4);
            objDtrNo_col = zeros(length(stim), 1);
            iCount   = 1;
            for iBlock = 1 : nBlock
                blc_i = find((subNo_col == (iSub - 1)) & (blockNo_col == (iBlock - 1)));
                trlNo_i   = trialNo_col(blc_i);
                tgtAng_i  = tgtAng_col(blc_i);
                dtrAng1_i = dtrAng1_col(blc_i);
                dtrAng2_i = dtrAng2_col(blc_i);
                dtrAng3_i = dtrAng3_col(blc_i);
                dtrNoCnt_i= dtrNoCnt_col(blc_i);

                load([subjDir, subjBv,  'mouseTraj_trials_in_Blc', num2str(iBlock), '.mat'], 'mouseTraj_trials');
                for jTrl = 1 : length(mouseTraj_trials)
                    trl_j     = find(trlNo_i == (jTrl - 1));
                    tgtAng_j  = tgtAng_i(trl_j);
                    dtrAng1_j = dtrAng1_i(trl_j);
                    dtrAng2_j = dtrAng2_i(trl_j);
                    dtrAng3_j = dtrAng3_i(trl_j);
                    objAngs   = [tgtAng_j, dtrAng1_j, dtrAng2_j, dtrAng3_j];
                    dtrNoCnt_j= dtrNoCnt_i(trl_j);
                    
                    timePass   = mouseTraj_trials{jTrl}(:, 4);
                    cueTp      = find(timePass <= 0.8);
                    cueTp_stay = cueTp(end); %% the final stay point before cue onset
                    objAngs_stay = objAngs(cueTp_stay, :);
                    objAng_col(iCount, :) = objAngs_stay;
                    objDtrNo_col(iCount)  = dtrNoCnt_j(1);
                    
                    objAngs_stay_left = objAngs_stay(1 : (dtrNoCnt_j + 1));
                    if ~any(isnan(objAngs_stay_left))
                        [~, minId] = min(abs(objAngs_stay_left));
                        if minId == 1
                            choiceId_ang(iCount) = 1; % 1: choose the target
                            
                        elseif minId ~= 1
                            choiceId_ang(iCount) = 0; % 0: choose the non-target
                            
                        end
                    end
                    iCount = iCount + 1;
                end
            end
            trl_counts = (1 : 1 : length(stim))'; 
            
            %% classify 'Yes' or 'No' according to angles in mouse trajectory or subjects' selections
            if isequal(fitWord, 'angleChoice')
                choiceId = choiceId_ang;
                if nanRemove == 1
                    trl_counts(isnan(choiceId)) = [];
                    stim(isnan(choiceId))       = [];
                    resp(isnan(choiceId), :)    = [];
                    objAng_col(isnan(choiceId), :)   = [];
                    objDtrNo_col(isnan(choiceId), :) = [];
                    %%% In the incorrect version, choiceId(isnan(choiceId)) =
                    %%% [] was run before these following 4 variables
                    delta_t_trials(isnan(choiceId), :)    = [];
                    assDist_mat(:, :, isnan(choiceId))    = [];
                    relation_history(:, :, isnan(choiceId)) = [];
                    posXY(isnan(choiceId), :, :)          = [];
                    respYes_trials_blc(isnan(choiceId))   = [];
                    mouseTraj_trials_blc(isnan(choiceId)) = [];
                    transStyle(isnan(choiceId))           = [];
                    rndOrHam(isnan(choiceId))             = [];
                    dtNum(isnan(choiceId))                = [];
                    choiceId(isnan(choiceId))             = [];
                end
                
            elseif isequal(fitWord, 'clickChoice')
                choiceId = choiceId_click;
                %%% !!!!!!!! remove the time-out trials !!!!!!!! %%%
                errTrials = find(errorIdx_Col == 1);
                errLeft   = zeros(length(errTrials), 1);
                for i = 1 : length(errTrials)
                    err_i     = errTrials(i);
                    respYes_i = respYes_trials_blc{err_i};
                    inCid     = find(respYes_i == -1);
                    if isempty(inCid) %% time-out trials
                        errLeft(i) = 1;
                    end
                end
                exTid = errTrials(errLeft == 1);
                
                trl_counts(exTid)  = [];
                stim(exTid)        = [];
                resp(exTid, :)     = [];
                posXY(exTid, :, :) = [];
                errorIdx_Col(exTid)= [];
                choiceId(exTid)    = [];
                respYes_trials_blc(exTid)   = [];
                mouseTraj_trials_blc(exTid) = [];
                delta_t_trials(exTid, :)    = [];
                assDist_mat(:, :, exTid)    = [];
                relation_history(:, :, exTid) = [];
                objAng_col(exTid, :)        = [];
                objDtrNo_col(exTid, :)      = [];

            end
            if isequal(trialWord, 'Hamiltonian')
                stim                 = stim(trl_counts > rndTrial);
                resp                 = resp(trl_counts > rndTrial, :);
                choiceId             = choiceId(trl_counts > rndTrial);
                delta_t_trials       = delta_t_trials(trl_counts > rndTrial, :);
                assDist_mat          = assDist_mat(:, :, trl_counts > rndTrial);
                relation_history     = relation_history(:, :, trl_counts > rndTrial);
                posXY                = posXY(trl_counts > rndTrial, :, :);
                respYes_trials_blc   = respYes_trials_blc(trl_counts > rndTrial);
                mouseTraj_trials_blc = mouseTraj_trials_blc(trl_counts > rndTrial);
                objAng_col           = objAng_col(trl_counts > rndTrial, :);
                objDtrNo_col         = objDtrNo_col(trl_counts > rndTrial, :);
            end
            transStyle_subj{iSub} = transStyle;
            rndOrHam_subj{iSub}   = rndOrHam;
            dtNum_subj{iSub}      = dtNum;

            %% fitting procedure
            nData_i = length(stim);
            if isequal(iteWord, 'naive')       %% learn the transition matrix from zero
                Aest_input = eye(states, states);
                nij_input  = zeros(states, states);
                M_input    = eye(states, states);
            end
            for iM = 1 : length(MList)
                Midx      = MList{iM};
                priorWord = priorList{iM}; %'InformPrior', 'nullPrior';

                if isequal(Midx, 'MaxEntropy')
                    transMat_input = Aest_input;
                else
                    transMat_input = M_input;
                end 
                if isequal(Midx, 'oneShotLinDecay') || isequal(Midx, 'oneShotExpDecay') ...
                        || isequal(Midx, 'SRoneShotLinDecay') || isequal(Midx, 'SRoneShotExpDecay') ...
                        || isequal(Midx, 'SRAssOneshotLin3Comp')
                    assDist_mat = relation_history;
                end

                [lliEst, paramsEst, lli_all] = ClusterRep_FitModels_choice(ExpWord, subID, stim, resp, states, choiceId, delta_t_trials, assDist_mat, fitWord, ...
                                               expMode_i, posXY, respYes_trials_blc, mouseTraj_trials_blc, objAng_col, objDtrNo_col, miniDmat, nFit, Midx, folder, refit, transMat_input, nij_input, iteWord, trialWord, priorWord);
                
                if isequal(Midx, 'SRheurProcessProb')
                    [~, M_trans, M_ass, p_choice_all, choiceID_all, lambda_trials] = SR_heurProcessProb_choice_gdChoice(paramsEst, stim, resp, states, choiceId, delta_t_trials, assDist_mat, expMode_i, posXY, respYes_trials_blc, ...
                                                            mouseTraj_trials_blc, objAng_col, objDtrNo_col, miniDmat, fitWord, transMat_input, nij_input, priorWord);
                elseif isequal(Midx, 'SRheurProcessTraceWin')
                    [~, M_trans, M_ass, p_choice_all, choiceID_all, lambda_trials] = SR_heurProcessTraceWin_choice_gdChoice(paramsEst, stim, resp, states, choiceId, delta_t_trials, assDist_mat, expMode_i, posXY, respYes_trials_blc, ...
                                                            mouseTraj_trials_blc, objAng_col, objDtrNo_col, miniDmat, fitWord, transMat_input, nij_input, priorWord);
                elseif isequal(Midx, 'SRheurProcessTrace')
                    [~, M_trans, M_ass, p_choice_all, choiceID_all, lambda_trials] = SR_heurProcessTrace_choice_gdChoice(paramsEst, stim, resp, states, choiceId, delta_t_trials, assDist_mat, expMode_i, posXY, respYes_trials_blc, ...
                                                            mouseTraj_trials_blc, objAng_col, objDtrNo_col, miniDmat, fitWord, transMat_input, nij_input, priorWord);
                elseif isequal(Midx, 'SRheurProcessTraceTrunc')
                    [~, M_trans, M_ass, p_choice_all, choiceID_all, lambda_trials] = SR_heurProcessTrace_choice_gdChoice_trunc(paramsEst, stim, resp, states, choiceId, delta_t_trials, assDist_mat, expMode_i, posXY, respYes_trials_blc, ...
                                                            mouseTraj_trials_blc, objAng_col, objDtrNo_col, miniDmat, fitWord, transMat_input, nij_input, priorWord);  
                elseif isequal(Midx, 'SRheurProcessTraceTruncAtt')
                    [~, M_trans, M_ass, p_choice_all, choiceID_all, lambda_trials] = SR_heurProcessTrace_choice_gdChoice_truncAtt(paramsEst, stim, resp, states, choiceId, delta_t_trials, assDist_mat, expMode_i, posXY, respYes_trials_blc, ...
                                                            mouseTraj_trials_blc, objAng_col, objDtrNo_col, miniDmat, fitWord, transMat_input, nij_input, priorWord);
                elseif isequal(Midx, 'SRheurProcessTraceTruncAttMod')
                    [~, M_trans, M_ass, p_choice_all, choiceID_all, lambda_trials] = SR_heurProcessTrace_choice_gdChoice_truncAttMod(paramsEst, stim, resp, states, choiceId, delta_t_trials, assDist_mat, expMode_i, posXY, respYes_trials_blc, ...
                                                            mouseTraj_trials_blc, objAng_col, objDtrNo_col, miniDmat, fitWord, transMat_input, nij_input, priorWord);
                end

                M_trans_subj(:, :, iSub, iM) = M_trans;
                M_ass_subj(:, :, iSub, iM)   = M_ass;
                lambda_trials_subj_nExp(iSub, trl_counts, iExp) = lambda_trials;

                %%% ------ quantifying lambda_trials: the weight after within and between transitions ------
                % ---- All trials ----
                lambda_mean_subj_nExp(iSub, 1, iExp) = nanmean(lambda_trials(transStyle == 0));
                lambda_mean_subj_nExp(iSub, 2, iExp) = nanmean(lambda_trials(transStyle == 1));
                % ---- The final 800 trials (mixture of random and hamiltonian walk) ----
                lambda_mean_subj_nExp(iSub, 3, iExp) = nanmean(lambda_trials(find(transStyle == 0 & trl_counts > rndTrial)));
                lambda_mean_subj_nExp(iSub, 4, iExp) = nanmean(lambda_trials(find(transStyle == 1 & trl_counts > rndTrial)));

                %%% ------ lambda_trials for transitions (within vs. between) and walks (random vs. hamiltonian) ------
                %%% only for the last 800 trials (mixture of random and hamiltonian walk)
                % ---- Random Walk ----
                lambda_mean_trans_walk_nExp(iSub, 1, 1, iExp) = nanmean(lambda_trials(find(transStyle == 0 & rndOrHam == 1 & trl_counts > rndTrial)));
                lambda_mean_trans_walk_nExp(iSub, 2, 1, iExp) = nanmean(lambda_trials(find(transStyle == 1 & rndOrHam == 1 & trl_counts > rndTrial)));
                % ---- Hamiltonian Walk ----
                lambda_mean_trans_walk_nExp(iSub, 1, 2, iExp) = nanmean(lambda_trials(find(transStyle == 0 & rndOrHam == 2 & trl_counts > rndTrial)));
                lambda_mean_trans_walk_nExp(iSub, 2, 2, iExp) = nanmean(lambda_trials(find(transStyle == 1 & rndOrHam == 2 & trl_counts > rndTrial)));
            end
            
            %%
            
        elseif isequal(expMode_i, 'key')
            % ......
 
        end
    end
end

%% behavioral paper, Figure xx: lambda_mean for all trials or the final 800 trials
dataFlg = 1;
if dataFlg == 1     % ------ all trials ------
    lambda_mean_plot = lambda_mean_subj_nExp(:, 1 : 2, :);
elseif dataFlg == 2 % ------ final 800 trials ------
    lambda_mean_plot = lambda_mean_subj_nExp(:, 3 : 4, :);
end

figKey = 1;  % 0: figure for presentation; 1: figure for AI.
if figKey == 0
    barLineWid = 2;
    errLineWid = 3;
    refLineWid = 1;
elseif figKey == 1
    barLineWid = 1;
    errLineWid = 1.5; %2;
    refLineWid = 0.5;
end
barPos = [1, 1.5];
for iExp = 1 : 3
    lambda_mean_plot_i = lambda_mean_plot(:, :, iExp);
    disp(['---------- ', expList{iExp}, ' ----------']);
    figure('Position', [100 100 200 120]), clf;

    for ib = 1 : 2 % within vs. between transitions
        [tAcc_avg, tAcc_sem] = Mean_and_Se(lambda_mean_plot_i, 1);
        plot(barPos, lambda_mean_plot_i, 'Color', [0.6, 0.6, 0.6], 'LineStyle', '-', 'LineWidth', 0.4); hold on;
        plot(barPos, tAcc_avg, 'Color', [0, 0, 0], 'LineStyle', '-', 'LineWidth', errLineWid); hold on;
        for iDm = 1 : 2
            errorbar(barPos(iDm), tAcc_avg(iDm), tAcc_sem(iDm), 'Color', 'k', 'LineStyle', 'none', 'LineWidth', errLineWid); hold on;
            plot(barPos(iDm), tAcc_avg(iDm), 'Marker', 'o', 'MarkerSize', 4.5, 'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', 'k', 'LineStyle', '-'); hold on;
        end
    end
    % ------ Statistical tests ------
    disp('======== \lambda~=0.5: within ========');
    [h, p, ci, stats] = ttest(lambda_mean_plot_i(:, 1), 0.5)
    disp('======== \lambda~=0.5: between ========');
    [h, p, ci, stats] = ttest(lambda_mean_plot_i(:, 2), 0.5)
    disp('======== \lambda: within vs. between ========');
    [h, p, ci, stats] = ttest(lambda_mean_plot_i(:, 1), lambda_mean_plot_i(:, 2))
    xlim([0.6, 1.9]);
    ylim([0, 1]);
    plot(xlim, [0.5, 0.5], 'k--', 'LineWidth', 0.6); hold on;
    if figKey == 0
        % ------For presentation------
        set(gca, 'LineWidth', 2);
        set(gca, 'FontSize', 15, 'FontWeight', 'bold', 'FontName', 'Arial');
        set(gca, 'XTick', '', 'XTickLabel', '');
        set(gca, 'YTick', 0 : 0.5 : 1, 'YTickLabel', 0 : 0.5 : 1);
    elseif figKey == 1
        % ------For Adobe Illustrator------
        set(gca, 'LineWidth', 0.8); % 0.8
        set(gca, 'FontSize', 10, 'FontWeight', 'bold', 'FontName', 'Arial');
        set(gca, 'XTick', [1, 1.5], 'XTickLabel', '');
        set(gca, 'YTick', 0 : 0.5 : 1, 'YTickLabel', {'', '', ''});
    end
    box off;
    ax = gca;
    save_name = ['Exp', num2str(iExp), '-', expList{iExp}, '-lambda-mean.png'];
    %exportgraphics(ax, save_name, 'Resolution', 600);

end

%% behavioral paper, Figure xx: lambda_mean for 2 transitions (within vs. between) times 2 walks (random vs. hamiltonian) or the final 800 trials
% lambda_mean_trans_walk_nExp = nan(subLen, 2, 2, length(expList));
figKey = 1;
if figKey == 0
    barLineWid = 2;
    errLineWid = 3;
    refLineWid = 1;
    indvLineW  = 1;
    markSize   = 6;
elseif figKey == 1
    barLineWid = 1;
    errLineWid = 2;
    refLineWid = 0.5;
    indvLineW  = 0.4;
    markSize   = 4.5;
end

barPos   = [1, 1.5; 1.7, 2.2];
for iExp = 1 : 3
    disp(['---------- ', expList{iExp}, ' ----------']);
    figure('Position', [100 100 260 120]), clf;

    lambda_trans_walk_iExp = lambda_mean_trans_walk_nExp(:, :, :, iExp);
    [accAvg_rOh, accSem_rOh] = Mean_and_Se(lambda_trans_walk_iExp, 1);
    accAvg_rOh = squeeze(accAvg_rOh);
    accSem_rOh = squeeze(accSem_rOh); % acc_trans_walk = zeros(subLen, 2, 2, length(expList));
    for i_rOh = 1 : 2 % Random vs. Hamiltonian Walk
        if i_rOh == 1
            walkWord = 'Radom';
        elseif i_rOh == 2
            walkWord = 'Hamiltonian';
        end
        barPos_i = barPos(i_rOh, :);
        %%% line plot
        plot(barPos_i, lambda_trans_walk_iExp(:, :, i_rOh), 'Color', [0.6, 0.6, 0.6], 'LineStyle', '-', 'LineWidth', indvLineW); hold on;
        plot(barPos_i, accAvg_rOh(:, i_rOh), 'Color', [0, 0, 0], 'LineStyle', '-', 'LineWidth', errLineWid); hold on;
        for i_tp = 1 : 2 % within vs. between-transition
            errorbar(barPos_i(i_tp), accAvg_rOh(i_tp, i_rOh), accSem_rOh(i_tp, i_rOh), 'Color', 'k', 'LineStyle', 'none', 'LineWidth', errLineWid); hold on;
            plot(barPos_i(i_tp), accAvg_rOh(i_tp, i_rOh), 'Marker', 'o', 'MarkerSize', markSize, 'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', 'k', 'LineStyle', '-'); hold on;
        end
        % ------ Statistical tests ------
        disp(['======== ', walkWord, '-within vs. between trans ========']);
        [h, p, ci, stats] = ttest(lambda_trans_walk_iExp(:, 1, i_rOh), lambda_trans_walk_iExp(:, 2, i_rOh))
    end
    xlim([0.6, 2.6]);
    ylim([0, 1]);
    plot(xlim, [0.5, 0.5], 'k--', 'LineWidth', 0.6); hold on;
    if figKey == 0
        % ------For presentation------
        set(gca, 'LineWidth', 2);
        set(gca, 'FontSize', 15, 'FontWeight', 'bold', 'FontName', 'Arial');
        set(gca, 'XTick', '', 'XTickLabel', '');
        set(gca, 'YTick', 0 : 0.5 : 1, 'YTickLabel', 0 : 0.5 : 1);
    elseif figKey == 1
        % ------For Adobe Illustrator------
        set(gca, 'LineWidth', 0.8);
        set(gca, 'FontSize', 10, 'FontWeight', 'bold', 'FontName', 'Arial');
        set(gca, 'XTick', [1, 1.5, 1.7, 2.2], 'XTickLabel', '');
        set(gca, 'YTick', 0 : 0.5 : 1, 'YTickLabel', '');
    end
    box off;
    ax = gca;
    save_name = ['Exp', num2str(iExp), '-', expList{iExp}, '-lambda-trans-walk.png'];
%     exportgraphics(ax, save_name, 'Resolution', 600);

end

%% Behavioral paper, Figure xx: weights trace across trials (bin every 100 trials)
% lambda_trials_subj_nExp = nan(subLen, nTrials, length(expList));
lambda_bin_subj_nExp = nan(subLen, BinL, length(expList));
for iExp = 1 : 3
    for iSub = 1 : size(lambda_trials_subj_nExp, 1)
        lambda_trials_iSub = squeeze(lambda_trials_subj_nExp(iSub, :, iExp));
        for iB = 1 : BinL
            trlIdx    = (iB - 1) * nBin + 1 : iB * nBin;
            lambda_iB = lambda_trials_iSub(trlIdx);
            lambda_bin_subj_nExp(iSub, iB, iExp) = nanmean(lambda_iB);
        end
    end
end

figKey = 1;
if figKey == 0
    barLineWid = 2;
    errLineWid = 3;
    refLineWid = 1;
    indvLineW  = 1;
    markSize   = 6;
elseif figKey == 1
    barLineWid = 1;
    errLineWid = 2;
    refLineWid = 0.5;
    indvLineW  = 0.4;
    markSize   = 6;
end

for iExp = 1 : length(expList)
    figure('Position', [100 100 260 120]), clf;
    hold on;
    
    lambdaAcc_ii = squeeze(lambda_bin_subj_nExp(:, :, iExp));
    [accAvg, accSem] = Mean_and_Se(lambdaAcc_ii, 1);
    errorbar(1 : 1 : BinL, accAvg, accSem, 'Color', [0.78, 0.50, 0.75], 'LineStyle', '-', 'LineWidth', errLineWid); hold on;

    for iB = 1 : BinL
        plot(iB, accAvg(iB), 'Marker', 'o', 'MarkerSize', markSize, 'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', [0.78, 0.50, 0.75], 'LineStyle', '-'); hold on;
    end
    xlim([0.5, BinL+0.5]);
    ylim([0, 1]);
    plot(xlim, [0.5, 0.5], 'k--', 'LineWidth', 0.6); hold on;
    if figKey == 0
        % ------For presentation------
        set(gca, 'LineWidth', 2);
        set(gca, 'FontSize', 15, 'FontWeight', 'bold', 'FontName', 'Arial');
        set(gca, 'XTick', 1 : 1 : BinL, 'XTickLabel', 1 : 1 : BinL);
        set(gca, 'YTick', 0 : 0.5 : 1, 'YTickLabel', 0 : 0.5 : 1);
    elseif figKey == 1
        % ------For Adobe Illustrator------
        set(gca, 'LineWidth', 0.8);
        set(gca, 'FontSize', 10, 'FontWeight', 'bold', 'FontName', 'Arial');
        set(gca, 'XTick', 1 : 1 : BinL, 'XTickLabel', '');
        set(gca, 'YTick', 0 : 0.5 : 1, 'YTickLabel', '');
    end
    box off;
    ax = gca;
    save_name = ['Exp', num2str(iExp), '-', expList{iExp}, '-lambda-bins.png'];
    exportgraphics(ax, save_name, 'Resolution', 600);
end

%% Behavioral paper, Supplementary Figure xx: weights trace across trials (all 1500 trials)
% lambda_trials_subj_nExp = nan(subLen, nTrials, length(expList));
figKey = 1;
if figKey == 0
    barLineWid = 2;
    errLineWid = 3;
    refLineWid = 1;
    indvLineW  = 1;
    markSize   = 6;
elseif figKey == 1
    barLineWid = 1;
    errLineWid = 1;
    refLineWid = 0.5;
    indvLineW  = 0.4;
    markSize   = 4.5;
end

for iExp = 1 : 3
    figure('Position', [100 100 600 200]), clf;
    [lambdaAvg, lambdaSe] = Mean_and_Se(lambda_trials_subj_nExp(:, :, iExp), 1);
    %errorbar(1 : 1 : 1500, lambdaAvg, lambdaSe, 'Color', 'r', 'LineStyle', '-', 'LineWidth', 1); hold on;
    shadedErrorBar(1 : 1 : 1500, lambdaAvg, lambdaSe, {'Color', [0.78, 0.50, 0.75], 'MarkerFaceColor', [0.78, 0.50, 0.75], 'LineStyle', '-', 'LineWidth', errLineWid}, 0.5); hold on;
    ylim([0, 1]);
    plot(xlim, [0.5, 0.5], 'k--', 'LineWidth', 0.6); hold on;
    if figKey == 0
        % ------For presentation------
        set(gca, 'LineWidth', 2);
        set(gca, 'FontSize', 15, 'FontWeight', 'bold', 'FontName', 'Arial');
        set(gca, 'XTick', '', 'XTickLabel', '');
        set(gca, 'YTick', 0 : 0.5 : 1, 'YTickLabel', 0 : 0.5 : 1);
    elseif figKey == 1
        % ------For Adobe Illustrator------
        set(gca, 'LineWidth', 0.8);
        set(gca, 'FontSize', 10, 'FontWeight', 'bold', 'FontName', 'Arial');
        set(gca, 'XTick', [1, 1.5, 1.7, 2.2], 'XTickLabel', '');
        set(gca, 'YTick', 0 : 0.5 : 1, 'YTickLabel', '');
    end
    box off;
    ax = gca;
    save_name = ['Exp', num2str(iExp), '-', expList{iExp}, '-lambda-trials.png'];
    exportgraphics(ax, save_name, 'Resolution', 600);

end


%%
addpath('tight_subplot/');

%% single subject
% lambda_trials_subj_nExp = nan(subLen, nTrials, length(expList));
for iExp = 1 : 3
    figure('Position', [100 100 1200 600]), clf;
    ha = tight_subplot(4, 6, [.02 .02], [.03 .01], [.02 .01]); %gap, marg_h, marg_w
    ha_i = 1;
    for iSub = 1 : subLen
        axes(ha(ha_i));
        plot(1 : 1 : 1500, lambda_trials_subj_nExp(iSub, :, iExp), 'Color', 'k', 'LineStyle', '-', 'LineWidth', 1); hold on;
        set(gca, 'LineWidth', 1);
        set(gca, 'FontSize', 8, 'FontWeight', 'bold', 'FontName', 'Arial');
        box off;
        ha_i = ha_i + 1;
    end
end








