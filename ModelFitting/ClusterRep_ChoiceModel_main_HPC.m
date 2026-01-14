function ClusterRep_ChoiceModel_main_HPC(iExp, iMode)
% -------- modified by XR @ Jan 13 2026 based on --------
% ClusterRep_ChoiceModel_main
% Run exactly one experiment (from expList) as one HPC job.
% Usage:
%   matlab -batch "ClusterRep_ChoiceModel_main_HPC(1,1)"
%   iExp  : 1..length(expList)
%   iMode : 1..length(expModeList)  (1=mouse, 2=key)

if nargin < 1 || isempty(iExp);  iExp  = 1; end
if nargin < 2 || isempty(iMode); iMode = 1; end

clearvars -except iExp iMode
clc

%%
addpath(genpath('/home/mpib/ren/rxj-neurocode/HierarchicalCluster/'));

%% Columns index for the clusterResult.mat
rndTrial = 700;
HamTrial = 800; 
nTrials  = rndTrial + HamTrial;

expModeList = {'mouse', 'key'};
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

%%
rng('shuffle');

%% =================== Model fitting procedure ===================
%% Experiment 1:
% I didn't regress out the confounding variables from RTs.
fitWord   = 'angleChoice'; % {'clickChoice', 'angleChoice'}
iteWord   = 'naive';     % {'naive', 'trained'}
trialWord = 'all'; % {'all', 'Hamiltonian'}
priorWord = 'InformPrior'; %'nullPrior';
states = 15;
nFit   = 100;
nSim   = 500;
expList = {'ImplicitExp', 'ExplicitExp', 'ImplicitRandExp'};
MList   = {'heurProcessTraceDecay', 'heurProcessTraceDecayParam'}; % 'objTP', 'MaxEntropy', 'MaxEntropy', 'SR', 'heurProcessProb', 'heurProcessTraceWin', 'heurProcessTrace', 'heurProcessTraceAvg', 'SRheurProcessProb', 'SRheurProcessTraceWin', 'SRheurProcessTrace', 'SRheurProcessTraceAvg', 'SRheurProcessTraceCluster', 'SRheurProcessTraceDistort'
folder  = '/home/mpib/ren/rxj-neurocode/HierarchicalCluster'; %pwd;
refit   = 1;
calFlg  = 0; % 0: one-time; 1: accumulated
nanRemove = 1; % 1: remove 'nan' from choiceId and other variales before importing them into fitting procedure

% ---- choose exactly one experiment ----
if iExp < 1 || iExp > numel(expList)
    error('iExp out of range. Must be 1..%d', numel(expList));
end
ExpWord = expList{iExp};
fprintf('=== Running ExpWord = %s (iExp=%d), expMode=%s (iMode=%d) ===\n', ExpWord, iExp, expMode_i, iMode);

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

%% fitting for single subject
for iSub = 1 : subLen
    %%
    subID  = subjLab{iSub};
    subjBv = subj_listBv{iSub};
    subjDir = [bhvDataDir, subjBv, '/'];
    %%
    if isequal(expMode_i, 'mouse')
        load([subjDir, subjBv, 'clusterResult_Blocks.mat'], 'clusterResult');
        load([subjDir, subjBv,  'respYes_trials_blc.mat'], 'respYes_trials_blc');
        load([subjDir, subjBv,  'mouseTraj_trials_blc.mat'], 'mouseTraj_trials_blc');
        %%% !!!!!!!! remove the time-out trials !!!!!!!! %%%
        stim = clusterResult(:, objTgt);
        resp = clusterResult(:, objResp);
        posXY        = zeros(length(stim), length(imgPosX), 2);
        posXY(:, :, 1) = clusterResult(:, imgPosX);
        posXY(:, :, 2) = clusterResult(:, imgPosY);
        errorIdx_Col   = clusterResult(:, errorIdx); %% incorrect click or time-out
        timeOIdx_Col   = clusterResult(:, timeOIdx);
        fast_Col       = clusterResult(:, fastIdx);
        exId           = find(errorIdx_Col ~= 0);
        choiceId_click = ones(length(stim), 1);
        choiceId_click(exId) = 0;
            
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
                cueTp      = find(timePass <= 0.8); %% (1) if response before cue onset, select the last time point; (2) if response after cue onset, select the time piont at 0.8s;
                cueTp_stay = cueTp(end); %% the final stay point before cue onset
                objAngs_stay = objAngs(cueTp_stay, :);
                objAng_col(iCount, :) = objAngs_stay;
                objDtrNo_col(iCount)  = dtrNoCnt_j(1);
                
                objAngs_stay_left = objAngs_stay(1 : (dtrNoCnt_j + 1));
                if ~any(isnan(objAngs_stay_left))
                    [~, minId] = min(abs(objAngs_stay_left));
                    if minId == 1
                        choiceId_ang(iCount) = 1;
                        
                    elseif minId ~= 1
                        choiceId_ang(iCount) = 0;
                        
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
                delta_t_trials(isnan(choiceId), :) = [];
                assDist_mat(:, :, isnan(choiceId)) = [];
                relation_history(:, :, isnan(choiceId)) = [];
                posXY(isnan(choiceId), :, :)          = [];
                respYes_trials_blc(isnan(choiceId))   = [];
                mouseTraj_trials_blc(isnan(choiceId)) = [];
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
        %% fitting procedure
        if isequal(iteWord, 'naive')       %% learn the transition matrix from zero
            Aest_input = eye(states, states);
            nij_input  = zeros(states, states);
            M_input    = eye(states, states);
        end
        %%
        for iM = 1 : length(MList)
            Midx = MList{iM};
            if isequal(Midx, 'MaxEntropy')
                transMat_input = Aest_input;
            else%if isequal(Midx, 'SR') || isequal(Midx, 'SRgd') || isequal(Midx, 'SRgdC') || isequal(Midx, 'SRgdClusterTrace') || isequal(Midx, 'SRgdClusterWin') || isequal(Midx, 'SRcw') || isequal(Midx, 'objTP')
                transMat_input = M_input;
            end
            if isequal(Midx, 'oneShotLinDecay') || isequal(Midx, 'oneShotExpDecay') ...
                    || isequal(Midx, 'SRoneShotLinDecay') || isequal(Midx, 'SRoneShotExpDecay') ...
                    || isequal(Midx, 'SRAssOneshotLin3Comp')
                assDist_mat = relation_history;
            end

            %%
            [lliEst, paramsEst, lli_all] = ClusterRep_FitModels_choice(ExpWord, subID, stim, resp, states, choiceId, delta_t_trials, assDist_mat, fitWord, ...
                expMode_i, posXY, respYes_trials_blc, mouseTraj_trials_blc, objAng_col, objDtrNo_col, miniDmat, nFit, Midx, folder, refit, transMat_input, nij_input, iteWord, trialWord, priorWord);
        end
        
        %%
        
    elseif isequal(expMode_i, 'key')
        load([subjDir, subjBv, 'clusterResult_BlocksEyeTrack.mat'], 'clusterResult');
        load([subjDir, subjBv,  'respYes_trials_blc.mat'], 'respYes_trials_blc');
        load([subjDir, subjBv,  'mouseTraj_trials_blc.mat'], 'mouseTraj_trials_blc');
        %clusterResult = clusterResult(701 : end, :);
        clusterResult(isnan(clusterResult(:, 5)), :) = []; %% specifically for subject 'nwq_5_f_23' in the 1st experiment
        
        stim         = clusterResult(:, objTgt);
        resp         = clusterResult(:, objResp);
        errorIdx_Col = clusterResult(:, errorIdx);
        fixF_Col     = clusterResult(:, fixFroze);
        exId  = find(errorIdx_Col ~= 0 | fixF_Col > 0);
        exKey = find(errorIdx_Col ~= 0 & fixF_Col == 0); %% only incorrect key or time-out
        exFix = find(fixF_Col > 0);                      %% fixation failure included
        
        %%% classify trials into 3 categories: (1) response after cue;
        %%% (2) incorrect response before cue; (3) correct response
        %%% before cue
        respRT_Col = clusterResult(:, respRT);
        cueOn_Col  = clusterResult(:, cueOn);
        trueChoice_trials = zeros(length(stim), 3);
        trueChoice_trials(respRT_Col > cueOn_Col, 1) = 1;
        trueChoice_trials((respRT_Col <= cueOn_Col) & (errorIdx_Col ~= 0), 2) = 1;
        trueChoice_trials((respRT_Col <= cueOn_Col) & (errorIdx_Col == 0), 3) = 1;
        if unique(sum(trueChoice_trials, 2)) ~= 1
            error('Something is wrong!');
        end
        choiceId = trueChoice_trials;

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
        trl_counts = (1 : 1 : length(stim))';

        %%% !!!!!!!! remove the time-out trials !!!!!!!! %%%
        %%% ambiguous choice %%%
        errTrials = find(errorIdx_Col == 1);
        errLeft   = zeros(length(errTrials), 1);
        for i = 1 : length(errTrials)
            err_i     = errTrials(i);
            respYes_i = respYes_trials_blc{err_i};
            inCid     = find(respYes_i == 3);
            if ~isempty(inCid) %% time-out trials
                errLeft(i) = 1;
            end
        end  
        exTid = errTrials(errLeft == 1);
        
        if isequal(trialWord, 'Hamiltonian')
            stim           = stim(trl_counts > rndTrial);
            resp           = resp(trl_counts > rndTrial, :);
            choiceId       = choiceId(trl_counts > rndTrial);
            delta_t_trials = delta_t_trials(trl_counts > rndTrial, :);
        end

        %% fitting procedure
        if isequal(iteWord, 'naive')       %% learn the transition matrix from zero
            Aest_input = eye(states, states);
            nij_input  = zeros(states, states);
            M_input    = eye(states, states);
            
        elseif isequal(iteWord, 'trained') %% continue to learn the transition matrix entailed from 1st behavioral experiment
            MaxEntropy_Dir = [folder, '/ClusterRep_ChoiceModel_Results/FixFailureAsIncorrect/', 'keypress', '/', 'MaxEntropy', '-ModelFits/'];
            SR_Dir         = [folder, '/ClusterRep_ChoiceModel_Results/FixFailureAsIncorrect/', 'keypress', '/', 'SR', '-ModelFits/'];
            
            load([MaxEntropy_Dir, subjBv, '-', 'MaxEntropy', '-', 'naive', '-pred.mat'], 'pCorr', 'Aest', 'nij');
            load([SR_Dir, subjBv, '-', 'SR', '-', 'naive', '-pred.mat'], 'pCorr', 'M');
            Aest_input = Aest;
            nij_input  = nij;
            M_input    = M;
            
        end
        %%
        for iM = 1 : length(MList)
            %%
            Midx = MList{iM};
            if isequal(Midx, 'MaxEntropy')
                transMat_input = Aest_input;
            elseif isequal(Midx, 'SR')
                transMat_input = M_input;
            elseif isequal(Midx, 'SRsymm')
                transMat_input = M_input;
            elseif isequal(Midx, 'heurProcess')
                transMat_input = nij_input;
            elseif isequal(Midx, 'objTP')
                transMat_input = rdm_modelRDM([], [], [], 'minimalDist');
            elseif isequal(Midx, 'heurTP')
                transMat_input = LynnNetwork();
            end
            %% 
            [lliEst, paramsEst, lli_all] = ClusterRep_FitModels_choice(ExpWord, subID, stim, resp, states, choiceId, delta_t_trials, [], ...
                                           expMode_i, [], [], [], [], [], [], nFit, Midx, folder, refit, transMat_input, nij_input, iteWord, trialWord, priorWord); 
        end
    end
end
end