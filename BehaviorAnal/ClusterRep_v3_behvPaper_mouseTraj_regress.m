% ClusterRep_v3_mouseTraj_regress
% write by rxj @ 09/02/2021
%
% Linear regression on angles at each time point and for each individual
% =====================================================================
% Revision by rxj @ 04/13/2022, the code includes following two contents:
% (1) Linear regression on RTs;
% (2) Logistic regression on binary choice at the final decision point or
% at each sampling time point

clear
clc

%%
addpath('tight_subplot/');
addpath(genpath('HierarchicalCluster/'));

% addpath(genpath('/home/ubuntu/Xiangjuan/MyExperiment/HierarchicalCluster/'));

%% Columns index for the clusterResult.mat
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


%% =================== Model fitting procedure ===================
%% Experiment 1:
circle_list = 0 : 1/30 : 1.5; %0 : 1/60 : 1.5;
states      = 15;
expList     = {'ImplicitExp', 'ExplicitExp', 'ImplicitRandExp'}; %
folder      = '/Users/cdlab/Nextcloud/MATLAB/MyExperiment/HierarchicalCluster'; %pwd;
%folder      = '/Volumes/My Passport/renxiangjuan/Nextcloud/MATLAB/MyExperiment/HierarchicalCluster';
%folder          ='/home/ubuntu/Xiangjuan/MyExperiment/HierarchicalCluster';
rndTrial = 700;
HamTrial = 800; 
nTrials  = rndTrial + HamTrial;
rndInter = 85;
HamInter = 15;
grayoffT   = 0.8;

tgtNum  = 1;
dtrNums = [1, 2, 3];
nodeNum = 15;
clsDef  = [1, 2, 3, 4, 5; ...  %% the 1st cluster
           6, 7, 8, 9, 10; ... %% the 2nd cluster
           11, 12, 13, 14, 15];%% the 3rd cluster
transMat = LynnNetwork();
miniDmat = modelRDM_minimalDist();
paramsName = {'trialCnt', 'trialType_1', 'fromNode_1', 'toNode_1', 'Dcue', 'dtrNoCnt', 'recencyL', 'fromNode_1:toNode_1'}; % 'trialMarkLast_1',  
corrParams = {'trialCnt', 'trialType', 'fromNode', 'toNode', 'Dcue', 'dtrNoCnt', 'recencyL', 'fromNode:toNode'}; % 'trialMarkLast', 
nParams = length(paramsName); % exclude intercept parameter
nCorrs  = length(corrParams);

%% Linear or logistic regression
% regressors:
% trialCnt     : trial number, log scale;
% trialType    : random or hamiltonian trial£»
% fromNode     : boundary node (1) or within node (0);
% toNode       : boundary node (1) or within node (0);
% dtrInCue     : one of the distractors belongs to the same cluster as cue (1);
% dtrInTgt     : one of the distractors belongs to the same cluster as target (1);
% Dcue         : the minimal distance between distractor and cue;
% Dtgt         : the minimal distance between distractor and target;
% dtrNoCnt     : distractor numbers;
% trialMarkLast: correct (1) or incorrect response in the last trial;
% recencyL     : recency;
colName_RT    = {'RTs', 'trialCnt', 'trialType', 'fromNode', 'toNode', 'Dcue', 'dtrNoCnt', 'recencyL'}; % ,'trialMarkLast'
colName_bc    = {'Choice', 'trialCnt', 'trialType', 'fromNode', 'toNode', 'Dcue', 'dtrNoCnt', 'recencyL'}; % binary choice
colName_angle = {'Angle', 'trialCnt', 'trialType', 'fromNode', 'toNode', 'Dcue', 'dtrNoCnt', 'recencyL'}; % absolute angle
subLen = 24;
betasRT_subj    = nan(nParams, subLen, length(expList));
corrParams_subj = nan(nCorrs, nCorrs, subLen, length(expList)); 
betasChoice_subj= nan(nParams, subLen, length(expList));
betasChoice_time_subj = nan(nParams, length(circle_list), subLen, length(expList));
betasAngle_time_subj  = nan(nParams, length(circle_list), subLen, length(expList));
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
        tInBlc = [230, 2*230, 2*230+240, rndTrial + (200 : 200 : HamTrial)];
        
    elseif isequal(ExpWord, 'ExplicitExp')
        bhvDataDir = [folder, '/FormalExp-ExplicitLearning-LynnNetwork-Results/'];
        subj_list  = {'hsp_1_m_21', 'hr_2_f_22', 'pxy_3_m_19', 'wyx_4_m_22', 'lml_5_f_19', 'lf_6_m_21', 'md_7_f_21', 'srz_8_m_20', 'fsq_9_f_18', 'caq_10_f_20', ...
                      'xjm_11_f_26', 'cjs_12_m_19', 'xxx_13_f_19', 'wsq_14_f_19', 'zzm_15_f_18', 'lsy_16_f_20', 'man_17_f_19', 'zxy_18_f_18', 'cyh_19_f_21', 'szn_20_f_19', ...
                      'sb_21_m_22', 'drq_22_f_18', 'cy_23_f_18', 'zhc_24_m_21'};
        subjLab    = {'hsp1', 'hr2', 'pxy3', 'wyx4', 'lml5', 'lf6', 'md7', 'srz8', 'fsq9', 'caq10', ...
                      'xjm11', 'cjs12', 'xxx13', 'wsq14', 'zzm15', 'lsy16', 'man17', 'zxy18', 'cyh19', 'szn20', ...
                      'sb21', 'drq22', 'cy23', 'zhc24'};
        tInBlc = [230, 2*230, 2*230+240, rndTrial + (200 : 200 : HamTrial)];
        
    elseif isequal(ExpWord, 'ImplicitRandExp')
        bhvDataDir = [folder, '/FormalExp-LynnNetwork-ImplicitRandom-Results/'];
        subj_list  = {'zyh_1_f_21', 'why_2_m_18', 'cr_3_f_19', 'zyx_4_f_25', 'wym_5_m_23', 'wd_6_f_18', 'lyh_7_f_21', 'zr_8_f_24', 'zyh_9_m_19', 'zzy_10_f_19', ...
                      'smq_11_f_25', 'sz_12_f_20', 'lzy_13_m_22', 'yxy_14_f_18', 'zxl_15_f_18', 'wqh_16_m_23', 'zxj_17_m_21', 'skx_18_f_20', 'zlh_19_m_24', 'gwt_20_f_23', ...
                      'lwn_21_f_18', 'lrp_22_m_21', 'sjj_23_f_18', 'xy_24_f_19'}; % , 'xr_18_f_23'
        subjLab    = {'zyh1', 'why2', 'cr3', 'zyx4', 'wym5', 'wd6', 'lyh7', 'zr8', 'lyh9', 'zzy10', ...
                      'smq11', 'sz12', 'lzy13', 'yxy14', 'zxl15', 'wqh16', 'zxj17', 'skx18', 'zlh19', 'gwt20', ...
                      'lwn21', 'lrp22', 'sjj23', 'xy24'}; % , 'xr18'
        tInBlc = [230, 2*230, 2*230+240, rndTrial + (200 : 200 : HamTrial)];
        
    end
    nBlock = length(tInBlc); %% seperate the total trials into 10 blocks
    trialsInBlc = zeros(nBlock, 2);
    trialsInBlc(:, 1) = [1, tInBlc(1 : end - 1) + 1]';
    trialsInBlc(:, 2) = tInBlc;
        
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
        trialMark_col = angleTraj_pd_subj.trialMark2;
    end
    subj_listBv = subj_list;
    subLen      = length(subj_list); 

    %% individual subject
    for SubIdx = 1 : subLen
        disp([ExpWord, '-subj', num2str(SubIdx)]);
        %%
        subjBv      = subj_listBv{SubIdx};
        subTit      = subjLab{SubIdx};
        subjDir     = [bhvDataDir, subjBv, '/'];
        load([subjDir, subjBv, 'clusterResult_Blocks.mat'], 'clusterResult');
        
        trials_Col   = (1 : 1 : nTrials)';
        respRT_Col   = clusterResult(:, respRT);
        errorIdx_Col = clusterResult(:, errorIdx);
        timeOIdx_Col = clusterResult(:, timeOIdx);
        bonus_Col    = clusterResult(:, bonusCn);
        rndHam_Col   = clusterResult(:, trialTp);
        fast_Col     = clusterResult(:, fastIdx);
        rndHam_Col(rndHam_Col == 1) = 0;
        rndHam_Col(rndHam_Col == 2) = 1;
        trialLab     = zeros(nTrials, 1);
        trialLab(errorIdx_Col ~= 1 | fast_Col ~= 1) = 1; %% feedback: correct or incorrect
        
        %%% within vs. between cluster transition RTs
        from_nodes = clusterResult(:, objTgt);
        to_nodes   = clusterResult(:, objResp(1));
        dt_nodes   = clusterResult(:, objResp);
        %%% whether one of the distractors belong to the same cluster as the cue (from node)
        lureIn_cue = zeros(nTrials, 1); % whether one of the distractors belong to the same cluster as the cue (from node)
        for ic = 1 : size(clsDef, 1)
            if ic == 1
                clsIdx = find(from_nodes >= 1 & from_nodes <= 5);
            elseif ic == 2
                clsIdx = find(from_nodes >= 6 & from_nodes <= 10);
            elseif ic == 3
                clsIdx = find(from_nodes >= 11 & from_nodes <= 15);
            end
            for i_cls = 1 : length(clsIdx)
                lureFind = arrayfun(@(x) ismember(x, clsDef(ic, :)), dt_nodes(clsIdx(i_cls), 2 : end));
                lureIn_cue(clsIdx(i_cls), 1) = length(find(lureFind == 1));
            end
        end
        %%% whether one of the distractors belong to the same cluster as the target (to node)
        lureIn_tgt = zeros(nTrials, 1); % whether one of the distractors belong to the same cluster as the target (to node)
        for ic = 1 : size(clsDef, 1)
            if ic == 1
                clsIdx = find(to_nodes >= 1 & to_nodes <= 5);
            elseif ic == 2
                clsIdx = find(to_nodes >= 6 & to_nodes <= 10);
            elseif ic == 3
                clsIdx = find(to_nodes >= 11 & to_nodes <= 15);
            end
            for i_cls = 1 : length(clsIdx)
                lureFind = arrayfun(@(x) ismember(x, clsDef(ic, :)), dt_nodes(clsIdx(i_cls), 2 : end));
                lureIn_tgt(clsIdx(i_cls), 1) = (length(find(lureFind == 1)) >= 1);
            end
        end
        %%% option numbers, including target
        dtNum = arrayfun(@(x) length(find(dt_nodes(x, :)~=0)), (1 : nTrials)');
        %%% fromNode category: within- or between-cluster
        fromNodeLab = zeros(nTrials, 1);
        fromNodeLab(from_nodes == 1 | from_nodes == 5 | from_nodes == 6 | from_nodes == 10 | from_nodes == 11 | from_nodes == 15) = 1;
        %%% toNode category: within- or between-cluster
        toNodeLab = zeros(nTrials, 1);
        toNodeLab(to_nodes == 1 | to_nodes == 5 | to_nodes == 6 | to_nodes == 10 | to_nodes == 11 | to_nodes == 15) = 1;
        
        %%% minimal distance: miniDmat
        minD_cue = zeros(nTrials, 1);
        minD_tgt = zeros(nTrials, 1);
        for iT = 1 : nTrials
            %%% distance between distractors and fromNode (cue)
            minD_cue_i   = miniDmat((dt_nodes(iT, 2 : dtNum(iT)))', from_nodes(iT));
            minD_cue(iT) = nanmean(minD_cue_i);
            %%% distance between distractors and toNode (target)
            minD_tgt_i   = miniDmat((dt_nodes(iT, 2 : dtNum(iT)))', to_nodes(iT));
            minD_tgt(iT) = nanmean(minD_tgt_i);
        end
        
        %%% recency
        recency_nodes = zeros(nTrials, 1);
        for iT = 1 : nTrials
            rec_iT = find(to_nodes(iT) == to_nodes(1 : (iT - 1)));
            if isempty(rec_iT)
                recency_nodes(iT) = nodeNum;
            elseif ~isempty(rec_iT) && length(rec_iT) >= 1
                recency_nodes(iT) = iT - rec_iT(end);
            end
        end
        
        %% integrate the regressors
        glm_trial = [respRT_Col, log(trials_Col), rndHam_Col, fromNodeLab, toNodeLab, minD_cue, dtNum-1, recency_nodes]; % , trialLab
        %glm_trial = [glm_trial(2 : end, 1 : end - 1), glm_trial(1 : end - 1, end)]; 
        glm_trial(:, [1, 2, 6, 7, 8]) = zscore(glm_trial(:, [1, 2, 6, 7, 8]), 0, 1);
        %glm_trial(:, [1, 2, 6, 8]) = zscore(glm_trial(:, [1, 2, 6, 8]), 0, 1);
        
        %% correlation between regressors
        corr_cols = [glm_trial(:, 2 : end), glm_trial(:, 4).*glm_trial(:, 5)];
        corrParams_subj(:, :, SubIdx, iExp) = corr(corr_cols); 
        
        %% convert matrix into table
        lm_tbl = array2table(glm_trial, 'VariableNames', colName_RT);
        lm_tbl.trialType = categorical(lm_tbl.trialType); % nominal(RT_tbl.fromBdNode);
        lm_tbl.fromNode  = categorical(lm_tbl.fromNode);
        lm_tbl.toNode    = categorical(lm_tbl.toNode);
        %lm_tbl.dtrNoCnt  = categorical(lm_tbl.dtrNoCnt);
        %lm_tbl.trialMarkLast = categorical(lm_tbl.trialMarkLast);
        
        %% ==========Linear regression on RTs==========
        % colName_RT = {'RTs', 'trialCnt', 'trialType', 'fromNode', 'toNode', 'Dcue', 'dtrNoCnt', 'recencyL', 'trialMarkLast'};
        fit = fitlm(lm_tbl, 'RTs ~ trialCnt + trialType + fromNode + toNode + Dcue + dtrNoCnt + recencyL + fromNode*toNode'); % + trialMarkLast 
        fitNames = fit.CoefficientNames;
        betas = nan(nParams, 1);
        for iPs = 1 : nParams
            params_iPs = paramsName{iPs};
            if any(strcmp(fitNames, params_iPs))
                betas(iPs) = fit.Coefficients{params_iPs, 1}; % 1 is the column of the Estimate
                % betas = fit.Coefficients.Estimate;
            end
        end
        betasRT_subj(:, SubIdx, iExp) = betas;
      
        %% ==========Linear regression on binary choice at final decision time point==========
        load([subjDir, subjBv, 'mouseTraj_trials_blc.mat'], 'mouseTraj_trials_blc');
        stim = clusterResult(:, objTgt);
        resp = clusterResult(:, objResp);
        
        choiceId  = nan(length(stim), length(circle_list));
        choiceAng = nan(length(stim), length(circle_list));
        choiceId_trials = nan(length(stim), 1);
        iCount   = 1;
        for iBlock = 1 : nBlock
            blc_i = find((subNo_col == (SubIdx - 1)) & (blockNo_col == (iBlock - 1)));
            trlNo_i     = trialNo_col(blc_i);
            tgtAng_i    = tgtAng_col(blc_i);
            dtrAng1_i   = dtrAng1_col(blc_i);
            dtrAng2_i   = dtrAng2_col(blc_i);
            dtrAng3_i   = dtrAng3_col(blc_i);
            dtrNoCnt_i  = dtrNoCnt_col(blc_i);
            trialMark_i = trialMark_col(blc_i);
            
            ttBlcIdx = (trialsInBlc(iBlock, 1) : trialsInBlc(iBlock, 2))';
            mouseTraj_trials = mouseTraj_trials_blc(ttBlcIdx);
            for jTrl = 1 : length(mouseTraj_trials)
                trl_j     = find(trlNo_i == (jTrl - 1));
                tgtAng_j  = tgtAng_i(trl_j);
                dtrAng1_j = dtrAng1_i(trl_j);
                dtrAng2_j = dtrAng2_i(trl_j);
                dtrAng3_j = dtrAng3_i(trl_j);
                objAngs   = [tgtAng_j, dtrAng1_j, dtrAng2_j, dtrAng3_j];
                dtrNoCnt_j  = dtrNoCnt_i(trl_j);
                trialMark_j = trialMark_i(trl_j);
                trialMark   = unique(trialMark_j);
                
                timePass  = mouseTraj_trials{jTrl}(:, 4);
                for iTp = 1 : length(circle_list)
                    time_i = circle_list(iTp);
                    cueTp  = find(timePass <= time_i);
                    if ~isempty(cueTp)
                        cueTp_stay = cueTp(end); %% the final stay point before cue onset
                        if cueTp_stay > length(objAngs) %% ????????
                            objAngs_stay = objAngs(end, :);
                        else
                            objAngs_stay = objAngs(cueTp_stay, :);
                        end
                        
                        objAngs_stay_left = objAngs_stay(1 : (unique(dtrNoCnt_j) + 1));
                        if ~any(isnan(objAngs_stay_left))
                            choiceAng(iCount, iTp) = objAngs_stay_left(1);
                            
                            [~, minId] = min(abs(objAngs_stay_left));
                            if minId == 1
                                choiceId(iCount, iTp) = 1;
                                
                            elseif minId ~= 1
                                choiceId(iCount, iTp) = 0;
                                
                            end
                        end
                    end
                end
                
                %%
                %%% the choice when subjects responded before cue or at the time point where cue began to appear
                if trialMark == 0     %% correct response before cue
                    choiceId_trials(iCount) = 1;
                elseif trialMark == 2 %% incorrect response before cue
                    choiceId_trials(iCount) = 0;
                else                          %% response after cue
                    cueTp = find(circle_list <= grayoffT); %find(timePass <= grayoffT);
                    if ~isempty(cueTp)
                        cueTp_stay = cueTp(end);
                    end
                    choiceId_trials(iCount) = choiceId(iCount, cueTp_stay);
                end
                iCount = iCount + 1;                
                %%
%                 timePass   = mouseTraj_trials{jTrl}(:, 4);
%                 cueTp      = find(timePass <= 0.8); %% (1) if response before cue onset, select the last time point; (2) if response after cue onset, select the time piont at 0.8s;
%                 cueTp_stay = cueTp(end); %% the final stay point before cue onset
%                 objAngs_stay = objAngs(cueTp_stay, :);
%                 
%                 objAngs_stay_left = objAngs_stay(1 : (unique(dtrNoCnt_j) + 1));
%                 if ~any(isnan(objAngs_stay_left))
%                     [~, minId] = min(abs(objAngs_stay_left));
%                     if minId == 1
%                         choiceId_trials(iCount) = 1;
%                         
%                     elseif minId ~= 1
%                         choiceId_trials(iCount) = 0;
%                         
%                     end
%                 end
%                 iCount = iCount + 1;
                
            end
        end
        
        %% fitting
        logit_trial = glm_trial; % remove trialMarkLast
        logit_trial(:, 1)   = choiceId_trials;
        logit_trial(isnan(logit_trial(:, 1)), :) = [];
        
        %% convert matrix into table
        lm_tbl = array2table(logit_trial, 'VariableNames', colName_bc);
        lm_tbl.trialType = categorical(lm_tbl.trialType); % nominal(RT_tbl.fromBdNode);
        lm_tbl.fromNode  = categorical(lm_tbl.fromNode);
        lm_tbl.toNode    = categorical(lm_tbl.toNode);
        %lm_tbl.dtrNoCnt  = categorical(lm_tbl.dtrNoCnt);
        %lm_tbl.trialMarkLast = categorical(lm_tbl.trialMarkLast);
        
        %% fitting
        modelspec = 'Choice ~ trialCnt + trialType + fromNode + toNode + Dcue + dtrNoCnt + recencyL + fromNode*toNode'; % trialMarkLast + 
        opts = statset('glmfit');
        opts.MaxIter = nParams*1000;
        fit = fitglm(lm_tbl, modelspec, 'Distribution', 'binomial', 'link', 'logit', 'options', opts);
        fitNames = fit.CoefficientNames;
        betas = nan(nParams, 1);
        for iPs = 1 : nParams
            params_iPs = paramsName{iPs};
            if any(strcmp(fitNames, params_iPs))
                betas(iPs) = fit.Coefficients{params_iPs, 1}; % 1 is the column of the Estimate
                % betas = fit.Coefficients.Estimate;
            end
        end
        betasChoice_subj(:, SubIdx, iExp) = betas;
        
        %% Logistic regression on binary choice at each sampling time point
        betasChoice_time = nan(nParams, length(circle_list));
        for iTime = 1 : length(circle_list)
            choiceId_i = choiceId(:, iTime);
            logit_trial = glm_trial;
            logit_trial(:, 1) = choiceId_i;
            logit_trial(isnan(logit_trial(:, 1)), :) = [];
            
            if size(logit_trial, 1) > (3 * nParams)
                %%% convert matrix into table
                lm_tbl = array2table(logit_trial, 'VariableNames', colName_bc);
                lm_tbl.trialType = categorical(lm_tbl.trialType); % nominal(RT_tbl.fromBdNode);
                lm_tbl.fromNode  = categorical(lm_tbl.fromNode);
                lm_tbl.toNode    = categorical(lm_tbl.toNode);
                %lm_tbl.dtrNoCnt  = categorical(lm_tbl.dtrNoCnt);
                %%% fitting
                modelspec = 'Choice ~ trialCnt + trialType + fromNode + toNode + Dcue + dtrNoCnt + recencyL + fromNode*toNode';
                opts = statset('glmfit');
                opts.MaxIter = nParams*1000;
                fit = fitglm(lm_tbl, modelspec, 'Distribution', 'binomial', 'link', 'logit', 'options', opts); % 
                fitNames = fit.CoefficientNames;
                betas = nan(nParams, 1);
                for iPs = 1 : nParams
                    params_iPs = paramsName{iPs};
                    if any(strcmp(fitNames, params_iPs))
                        betas(iPs) = fit.Coefficients{params_iPs, 1}; % 1 is the column of the Estimate
                        % betas = fit.Coefficients.Estimate;
                    end
                end
                betasChoice_time(:, iTime) = betas;
            end
        end
        betasChoice_time_subj(:, :, SubIdx, iExp) = betasChoice_time; 
        
        %% Linear regression on absolute angle of target at each sampling time point
        betaAngle_time = nan(nParams, length(circle_list));
        for iTime = 1 : length(circle_list)
            choiceId_i = abs(choiceAng(:, iTime));
            logit_trial = glm_trial;
            logit_trial(:, 1) = choiceId_i;
            logit_trial(isnan(logit_trial(:, 1)), :) = [];
            logit_trial(:, 1) = zscore(logit_trial(:, 1));
            
            if size(logit_trial, 1) > (3 * nParams)
                %%% convert matrix into table
                lm_tbl = array2table(logit_trial, 'VariableNames', colName_angle);
                lm_tbl.trialType = categorical(lm_tbl.trialType); % nominal(RT_tbl.fromBdNode);
                lm_tbl.fromNode  = categorical(lm_tbl.fromNode);
                lm_tbl.toNode    = categorical(lm_tbl.toNode);
                %%% fitting
                fit = fitlm(lm_tbl, 'Angle ~ trialCnt + trialType + fromNode + toNode + Dcue + dtrNoCnt + recencyL + fromNode*toNode'); % + trialMarkLast
                fitNames = fit.CoefficientNames;
                betas = nan(nParams, 1);
                for iPs = 1 : nParams
                    params_iPs = paramsName{iPs};
                    if any(strcmp(fitNames, params_iPs))
                        betas(iPs) = fit.Coefficients{params_iPs, 1}; % 1 is the column of the Estimate
                        % betas = fit.Coefficients.Estimate;
                    end
                end
                betaAngle_time(:, iTime) = betas;
            end
        end
        betasAngle_time_subj(:, :, SubIdx, iExp) = betaAngle_time;
        
    end
    % trialType_db_i(strcmp(trialType_i, 'random')) = 1;
    % save([angleDir, 'betas_subj_', dataFrom, '.mat'], 'betas_subj');

end
%% color settings
colorSets = [0.98, 0.72, 0.69; ...
            0.97, 0.85, 0.67; ...
            0.33, 0.73, 0.83; ...
            0.72, 0.80, 0.88; ...
            0.82, 0.92, 0.78; ...
            0.86, 0.80, 0.89; ...
            0.78, 0.50, 0.75; ...
            0.75, 0.56, 0; ...
            0.40, 0.40, 0.40; ...
            0.54, 0.67, 0.20];
colorSets_d = colorSets;
colorSets_d(2, :) = [];

%% corrrelation of the regressors
%%% one-sample t test for each correlation coefficient: corrParams_subj = nan(nParams, nParams, subLen, length(expList)); 
corrSig_exp    = ones(nCorrs, nCorrs, length(expList));
for iExp = 1 : length(expList)
    for i = 1 : (nCorrs - 1)
        for j = (i + 1) : nCorrs
            [h, p, ci, stats] = ttest(squeeze(corrParams_subj(j, i, :, iExp)), 0, 'tail', 'both');
            corrSig_exp(j, i, iExp) = p;
        end
    end
end
%%% heatmap for correlation and p-value
for iExp = 1 : length(expList)
    corr_iExp = corrParams_subj(:, :, :, iExp);
    corr_iExp_avg = nanmean(corr_iExp, 3);
    corr_iExp_avg(1 : (nCorrs + 1) : end) = 0;
    %%% correlations
    figure('Position', [100 100 350 300]), clf;
    clims = [-1, 1];
    imagesc([1, nCorrs], [1, nCorrs], corr_iExp_avg, clims); 
    %colormap(gray(256));
    [colorCohr]=cbrewer('div', 'RdBu', 100, 'cubic');
    colormap(colorCohr(end : -1 : 1, :));
    colorbar;
    axis xy;
    set(gca, 'FontSize', 22, 'FontWeight', 'bold', 'FontName', 'Arial');
    %set(gca, 'XTick', 1 : 1 : nParams, 'XTickLabel', '', 'YTick', 1 : 1 : nParams, 'YTickLabel', '');
    
    %%% labeling the significance
    for i = 1 : (nCorrs - 1)
        for j = (i + 1) : nCorrs
            pValue = corrSig_exp(j, i, iExp);
            if pValue < 0.001
                text(j, i, '***', 'HorizontalAlignment', 'Center', 'FontName', 'Arial', 'FontSize', 10);
            elseif pValue < 0.01 && pValue >=0.001
                text(j, i, '**', 'HorizontalAlignment', 'Center', 'FontName', 'Arial', 'FontSize', 10);
            elseif pValue < 0.05 && pValue >=0.01
                text(j, i, '*', 'HorizontalAlignment', 'Center', 'FontName', 'Arial', 'FontSize', 10);
            end
        end
    end
end

%% betas estimated from RTs or choices
dataFlg = 1; %input('Which data set, 0-RTs, 1-choice: ');
if dataFlg == 0
    betas_data = betasRT_subj;
elseif dataFlg == 1
    betas_data = betasChoice_subj;
    %betas_data(11, :, [1, 3]) = nan;
end

%% statistical testing of the betas
betas_pval_exp = nan(nParams, length(expList));
for iExp = 1 : length(expList)
    for iParams = 1 : nParams
        betas_ii = squeeze(betas_data(iParams, :, iExp));
        [h, p, ci, stats] = ttest(betas_ii, 0, 'tail', 'both');
        betas_pval_exp(iParams, iExp) = p;
    end
end

%% plotting the betas for each parameters: as barplot
%paramsName = {'trialCnt', 'trialType_1', 'fromNode_1', 'toNode_1', 'dtrInCue_1', 'dtrInTgt_1', 'Dcue', 'Dtgt', 'dtrNoCnt', 'recencyL', 'trialMarkLast_1', 'fromNode_1:toNode_1', 'dtrInCue_1:dtrInTgt_1'};
%betasRT_subj     = nan(nParams, subLen, length(expList));
%betasChoice_subj = nan(nParams, subLen, length(expList));
for iExp = 1 : length(expList)
    betas_i = betas_data(:, :, iExp);
    [bAvg, bSem] = Mean_and_Se(betas_i, 2);
    %bAvg = median(betas_i, 2, 'omitnan');
    pval_i = betas_pval_exp(:, iExp);
    [~, ~, ~, pval_i] = fdr_bh(pval_i, 0.05, 'pdep'); % method: 'dep', 'pdep'
    nParams = length(bAvg);
    
    figure('Position', [100 100 550 200]), clf;
    for iParams = 1 : nParams
        b = bar(iParams, bAvg(iParams)); hold on;
        b.FaceColor = [1, 1, 1];
        if pval_i(iParams) < 0.001
            b.EdgeColor = [1, 0, 0];
        elseif pval_i(iParams) >= 0.001 && pval_i(iParams) < 0.01
            b.EdgeColor = colorSets(1, :);
        elseif pval_i(iParams) >= 0.01 && pval_i(iParams) < 0.05
            b.EdgeColor = colorSets(2, :);
        else
            b.EdgeColor = [0, 0, 0];
        end
        b.LineWidth = 3;
        errorbar(iParams, bAvg(iParams), bSem(iParams), 'Color', 'k', 'LineStyle', '-', 'LineWidth', 4); hold on;
        plot(iParams, bAvg(iParams), 'Marker', '.', 'MarkerSize', 20, 'Color', 'k', 'LineStyle', 'none'); hold on;
    end
    set(gca, 'FontSize', 20, 'FontWeight', 'Bold', 'LineWidth', 2);
    set(gca, 'XTick', 1 : 1 : nParams, 'XTickLabel', '');
    xlim([0, nParams+1]);
    if dataFlg == 0
        set(gca, 'YTick', -0.4 : 0.2 : 0.5, 'YTickLabel', -0.4 : 0.2 : 0.5);
        ylim([-0.4, 0.5]);
    elseif dataFlg == 1
        set(gca, 'YTick', -1 : 0.5 : 0.5, 'YTickLabel', -1 : 0.5 : 0.5);
        ylim([-1, 0.5]);
    end
    
    box off;
end


%% plotting the betas at each time point within a trial
% write by rxj @ 06/15/2022
timeFlg = input('which time betas, 0-choice, 1-angle: ');
if timeFlg == 0
    betas_time_subj = betasChoice_time_subj;
elseif timeFlg == 1
    betas_time_subj = betasAngle_time_subj;
end
%% statistical test for each variable at each time point
% betasChoice_time_subj = nan(nParams, length(circle_list), subLen, length(expList));
statMat_times_exp = nan(nParams, length(circle_list), length(expList));
for iExp = 1 : length(expList)
    for iP = 1 : nParams
        for iTime = 1 : length(circle_list)
            betas_i = squeeze(betas_time_subj(iP, iTime, :, iExp));
            [h, p, ci, stats] = ttest(betas_i, 0, 'tail', 'both');
            statMat_times_exp(iP, iTime, iExp) = p;
        end
    end
end

%% plotting 
for iExp = 1 : length(expList)
    figure('Position', [100 100 400 220]), clf;
    hold on;
    %%
    markPos = 0.99;
    for iP = 1 : nParams
        betas_iExp = squeeze(betas_time_subj(iP, :, :, iExp));  % length(circle_list) * subLen
        [coef_i, sem_i] = Mean_and_Se(betas_iExp, 2);
        pval_i = squeeze(statMat_times_exp(iP, :, iExp));
        
        %shadedErrorBar(circle_list, coef_i, sem_i, {'Color', colorSets(iP, :), 'MarkerFaceColor', colorSets(iP, :), 'LineStyle', '-', 'LineWidth', 2}, 0.5); hold on;
        plot(circle_list, coef_i, 'Color', colorSets(iP, :), 'MarkerFaceColor', colorSets(iP, :), 'LineStyle', '-', 'LineWidth', 2); hold on;
        [~, ~, ~, adj_p] = fdr_bh(pval_i, 0.05, 'pdep'); % method: 'dep', 'pdep'
        
        for iCir = 1 : length(circle_list)
            pval_j = adj_p(iCir);
            if pval_j < 0.01
                plot(circle_list(iCir), markPos, 'Marker', '.', 'MarkerSize', 13, 'Color', colorSets(iP, :), 'MarkerFaceColor', [1,1,1], 'LineStyle', 'none', 'LineWidth', 1.5); hold on;
            elseif pval_j < 0.05 && pval_j >= 0.01
                plot(circle_list(iCir), markPos, 'Marker', 'diamond', 'MarkerSize', 4, 'Color', colorSets(iP, :), 'LineWidth', 1.5); hold on;
            end
        end
        markPos = markPos - 0.05;
    end
    xlim([0.1, 1.2]);
    plot(xlim, [0, 0], 'k--', 'LineWidth', 1); hold on;
    if iExp == 1 || iExp == 3
        plot([0.8, 0.8], ylim, 'k:', 'LineWidth', 1); hold on;
    end
    if timeFlg == 0
        ylim([-1, 1]);
    elseif timeFlg == 1
        ylim([-0.2, 0.2]);
    end
    axis xy;
    set(gca, 'LineWidth', 2);
    set(gca, 'FontSize', 18, 'FontWeight', 'bold', 'FontName', 'Arial');
    box off;
end










%% Old codes
%% statistical testing between learning and testing stage
stageWord = {'learning', 'testing'};
lme_coef_stage_exp   = zeros(nParams, length(circle_list), length(stageWord), length(expList));
lme_coefSe_stage_exp = zeros(nParams, length(circle_list), length(stageWord), length(expList));
lme_pval_stage_exp   = zeros(nParams, length(circle_list), length(expList));
for iExp = 1 : length(expList)
    ExpWord = expList{iExp};
    if isequal(ExpWord, 'ImplicitExp')
        bhvDataDir = [folder, '/FormalExp-LynnNetwork-Results/']; 
    elseif isequal(ExpWord, 'ExplicitExp')
        bhvDataDir = [folder, '/FormalExp-ExplicitLearning-LynnNetwork-Results/'];
    elseif isequal(ExpWord, 'ImplicitRandExp')
        bhvDataDir = [folder, '/FormalExp-LynnNetwork-ImplicitRandom-Results/'];
    end
    angleDir = [bhvDataDir, 'angleCal/'];
    for iStg = 1 : length(stageWord)
        dataFrom = stageWord{iStg};
        load([angleDir, 'betas_subj_', dataFrom, '.mat'], 'betas_subj');
        if iStg == 1
            betas_subj_learn = betas_subj;
        elseif iStg == 2
            betas_subj_test  = betas_subj;
        end
        [betasAvg, betasSem] = Mean_and_Se(betas_subj, 1);
        lme_coef_stage_exp(:, :, iStg, iExp)   = (squeeze(betasAvg))'; % betas_subj = nan(subLen, length(circle_list), nParams);
        lme_coefSe_stage_exp(:, :, iStg, iExp) = (squeeze(betasSem))';
    end

    for iParams = 1 : nParams
        for iTime = 1 : length(circle_list)
            betas_ii = betas_subj_learn(:, iTime, iParams);
            betas_jj = betas_subj_test(:, iTime, iParams);
            [h, p, ci, stats] = ttest(betas_ii, betas_jj);
            lme_pval_stage_exp(iParams, iTime, iExp) = p;
        end
    end
end

%% plotting each beta time series in different stages: learning and testing
xLoc = 0.95;
for i = 1 : length(expList)
    coef_i = lme_coef_stage_exp(:, :, :, i);
    sem_i  = lme_coefSe_stage_exp(:, :, :, i);
    pval_i = lme_pval_stage_exp(:, :, i);
    
    figure('Position', [100 100 800 450]), clf;
    ha = tight_subplot(3, 3, [.05 .05], [.1 .02], [.05 .05]);
    ha_i = 1;
    
    for j = 1 : size(coef_i, 1)
        axes(ha(ha_i));
        
        shadedErrorBar(circle_list, coef_i(j, :, 1), sem_i(j, :, 1), {'Color', colorUsed(j, :), 'MarkerFaceColor', colorUsed(j, :), 'LineStyle', '-', 'LineWidth', 2}, 0.5); hold on;
        shadedErrorBar(circle_list, coef_i(j, :, 2), sem_i(j, :, 2), {'Color', colorUsed(j, :), 'MarkerFaceColor', colorUsed(j, :), 'LineStyle', ':', 'LineWidth', 2}, 0.5); hold on;
        %plot(circle_list, coef_i(j, :, 1), 'Color', colorUsed(j, :), 'LineStyle', '-', 'LineWidth', 2); hold on;
        %plot(circle_list, coef_i(j, :, 2), 'Color', colorUsed(j, :), 'LineStyle', ':', 'LineWidth', 2); hold on;
        [~, ~, ~, adj_p] = fdr_bh((pval_i(j, :))', 0.05, 'pdep'); % method: 'dep', 'pdep'
        ylimit = ylim;
        for iCir = 1 : length(circle_list)
            pval_j = pval_i(j, iCir); %adj_p(iCir);
            if pval_j < 0.01
                plot(circle_list(iCir), xLoc * ylimit(end), 'Marker', '.', 'MarkerSize', 8, 'Color', colorUsed(j, :), 'MarkerFaceColor', colorUsed(j, :), 'LineStyle', 'none', 'LineWidth', 1.5); hold on;
            elseif pval_j < 0.05 && pval_j >= 0.01
                plot(circle_list(iCir), xLoc * ylimit(end), 'Marker', '>', 'MarkerSize', 4, 'Color', colorUsed(j, :), 'MarkerFaceColor', colorUsed(j, :), 'LineStyle', 'none', 'LineWidth', 1.5); hold on;
            end
            
        end
        %ylim([-0.2, 0.5]);
        plot(xlim, [0, 0], 'k--', 'LineWidth', 1); hold on;
        if i == 1 || i == 3
            plot([0.8, 0.8], ylim, 'k:', 'LineWidth', 1); hold on;
        end
        axis xy;
        set(gca, 'LineWidth', 2);
        set(gca, 'FontSize', 15, 'FontWeight', 'bold', 'FontName', 'Arial');
        box off;
        ha_i = ha_i + 1;
    end
end
















