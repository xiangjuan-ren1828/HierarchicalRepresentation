% ClusterRep_v3_behvPaper_statTest.m
% Modified by XR @ Jan 20 2026 based on the followings
% =========================================================================
% ClusterRep_v3_statTest
% write by rxj @ 11/03/2020
% =========================================================================
% ClusterRep_v2_statTest
% write by rxj @ 09/27/2020
% Codes for statistical test

%%
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
tgtNum  = 1;
dtrNums = [1, 2, 3];
nodeNum = 15;
clsDef  = [1, 2, 3, 4, 5; ...  %% the 1st cluster
           6, 7, 8, 9, 10; ... %% the 2nd cluster
           11, 12, 13, 14, 15];%% the 3rd cluster
transMat = LynnNetwork();
miniDmat = modelRDM_minimalDist();

%% 3 Experiments:
% I didn't regress out the confounding variables from RTs.
circle_list = 0 : 1/60 : 1.5;
states      = 15;
nFit        = 100;
expList     = {'ImplicitExp', 'ExplicitExp', 'ImplicitRandExp'}; 
folder      = '/Users/ren/Projects-NeuroCode/MyExperiment/HierarchicalCluster';

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

    % ------ read the trajectory-based data for all participants together ------
    if isequal(expMode_i, 'mouse') %% read the mouseTraj.csv
        angleDir = [bhvDataDir, 'angleCal/'];
        angleTraj_pd_subj = readtable([angleDir, 'angleTraj_pd_subj.csv']);
        %%% variables
        subNo_col    = angleTraj_pd_subj.subNo;
        trialCnt_col = angleTraj_pd_subj.trialCnt;
        blockNo_col  = angleTraj_pd_subj.blockNo;
        trialNo_col  = angleTraj_pd_subj.trialNo;
        tgtAng_col   = angleTraj_pd_subj.tgtAng;
        dtrAng1_col  = angleTraj_pd_subj.dtrAng1;
        dtrAng2_col  = angleTraj_pd_subj.dtrAng2;
        dtrAng3_col  = angleTraj_pd_subj.dtrAng3;
        dtrNoCnt_col = angleTraj_pd_subj.dtrNoCnt; %% 1, 2, 3 
    end
    subj_listBv = subj_list;
    subLen      = length(subj_list);

    %% ----- Loop over participants to get the corresponding data for the subsequent GLMM/LMM analysis ----- 

    glmm_subj = [];
    for SubIdx = 1 : subLen
        disp(['======', ExpWord, '-subj', num2str(SubIdx), '======']);
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
        %%% within vs. between cluster transtion
        from_nodes = clusterResult(:, objTgt);
        to_nodes   = clusterResult(:, objResp(1));
        dt_nodes   = clusterResult(:, objResp);
        transIn    = zeros(nTrials, 3);
        lureIn     = zeros(nTrials, 1); % whether there are two items selected from the same cluster
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
            for i_cls = 1 : length(clsIdx)
                lureFind = arrayfun(@(x) ismember(x, clsDef(ic, :)), dt_nodes(clsIdx(i_cls), 2 : end));
                lureIn(clsIdx(i_cls), 1) = length(find(lureFind == 1));
            end
        end
        transStyle = (transIn(:, 1) == transIn(:, 2));

        %%% mouse trajectory trials
        load([subjDir, subjBv,  'respYes_trials_blc.mat'], 'respYes_trials_blc');
        load([subjDir, subjBv,  'mouseTraj_trials_blc.mat'], 'mouseTraj_trials_blc');

        %%% choice accuracy
        stim = clusterResult(:, objTgt);
        resp = clusterResult(:, objResp);
        %%% extract the angles
        choiceId_ang  = nan(length(stim), 1);
        choiceId      = nan(length(stim), length(circle_list));
        dtrNoCnt_iSub = nan(length(stim), 1);
        iCount = 1;
        for iBlock = 1 : nBlock
            blc_i = find((subNo_col == (SubIdx - 1)) & (blockNo_col == (iBlock - 1)));
            trlNo_i   = trialNo_col(blc_i);
            tgtAng_i  = tgtAng_col(blc_i);
            dtrAng1_i = dtrAng1_col(blc_i);
            dtrAng2_i = dtrAng2_col(blc_i);
            dtrAng3_i = dtrAng3_col(blc_i);
            dtrNoCnt_i= dtrNoCnt_col(blc_i);

            ttBlcIdx = (trialsInBlc(iBlock, 1) : trialsInBlc(iBlock, 2))';
            mouseTraj_trials = mouseTraj_trials_blc(ttBlcIdx);

            for jTrl = 1 : length(mouseTraj_trials)
                trl_j     = find(trlNo_i == (jTrl - 1));
                tgtAng_j  = tgtAng_i(trl_j);
                dtrAng1_j = dtrAng1_i(trl_j);
                dtrAng2_j = dtrAng2_i(trl_j);
                dtrAng3_j = dtrAng3_i(trl_j);
                objAngs   = [tgtAng_j, dtrAng1_j, dtrAng2_j, dtrAng3_j];
                dtrNoCnt_j= dtrNoCnt_i(trl_j);
                dtrNoCnt_iSub(iCount) = unique(dtrNoCnt_j);

                timePass  = mouseTraj_trials{jTrl}(:, 4);
                %%% click time point (response before cue) or time point
                %%% before cue onset (response after cue)
                cueTp        = find(timePass <= 0.8);
                cueTp_stay   = cueTp(end); %% the final stay point before cue onset
                objAngs_stay = objAngs(cueTp_stay, :);
                objAngs_stay_left = objAngs_stay(1 : (unique(dtrNoCnt_j) + 1));
                if ~any(isnan(objAngs_stay_left))
                    [~, minId] = min(abs(objAngs_stay_left));
                    if minId == 1
                        choiceId_ang(iCount) = 1; % 1: choose the target

                    elseif minId ~= 1
                        choiceId_ang(iCount) = 0; % 0: choose the non-target
                    end
                end

                %%% every time point
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
                            [~, minId] = min(abs(objAngs_stay_left));
                            if minId == 1
                                choiceId(iCount, iTp) = 1;

                            elseif minId ~= 1
                                choiceId(iCount, iTp) = 0;

                            end
                        end
                    end
                end
                iCount = iCount + 1;
            end
        end

        %%% recency quantification
        recency_nodes = zeros(nTrials, 3);
        for iT = 1 : nTrials
            rec_iT = find(from_nodes(iT) == from_nodes(1 : (iT - 1)));
            if isempty(rec_iT)
                recency_nodes(iT, 1) = nodeNum;
                recency_nodes(iT, 2) = nodeNum;
                recency_nodes(iT, 3) = nodeNum;

            elseif ~isempty(rec_iT) && length(rec_iT) == 1
                recency_nodes(iT, 1) = iT - rec_iT(end);
                recency_nodes(iT, 2) = nodeNum;
                recency_nodes(iT, 3) = nodeNum;

            elseif ~isempty(rec_iT) && length(rec_iT) == 2
                recency_nodes(iT, 1) = iT - rec_iT(end);
                recency_nodes(iT, 2) = iT - rec_iT(end - 1);
                recency_nodes(iT, 3) = nodeNum;

            elseif ~isempty(rec_iT) && length(rec_iT) > 2
                recency_nodes(iT, 1) = iT - rec_iT(end);
                recency_nodes(iT, 2) = iT - rec_iT(end - 1);
                recency_nodes(iT, 3) = iT - rec_iT(end - 2);
            end

        end

        %%% Distractor numbers in each trial
        dtNum = arrayfun(@(x) length(find(dt_nodes(x, :)~=0)), (1 : nTrials)');

        %%% The objective distance between Cue and Target is always 1
        %%% 1) Cue-distractor distance; 2) Target-distractor distance
        cue_dtr_distance = nan(nTrials, 3); % 1/2/3 distractors
        tgt_dtr_distance = nan(nTrials, 3);
        for iT = 1 : nTrials
            dtr = resp(iT, 2 : dtNum(iT));          % distractors for this trial
            cue_dtr_distance(iT, 1 : numel(dtr)) = miniDmat(stim(iT), dtr);
            tgt_dtr_distance(iT, 1 : numel(dtr)) = miniDmat(resp(iT, 1), dtr);
        end

        %%% Whehter the cue, target and distractor were boundary nodes
        boundary_node_id = [1 5 6 10 11 15];
        cue_BdNode = zeros(nTrials, 1); %% whether the "from" (or cue) nodes are boundary nodes
        cue_BdNode(find(from_nodes == 1 | from_nodes == 5 | from_nodes == 6 | from_nodes == 10 | from_nodes == 11 | from_nodes == 15)) = 1;
        tgt_BdNode = zeros(nTrials, 1); %% whether the "to" (or target) nodes are boundary nodes
        tgt_BdNode(find(to_nodes == 1 | to_nodes == 5 | to_nodes == 6 | to_nodes == 10 | to_nodes == 11 | to_nodes == 15)) = 1;
        dtr_BdNode = zeros(nTrials, 3); %% whether the distractor nodes are boundary nodes
        dt_nodes_only = dt_nodes(:, 2 : end);
        mask = ismember(dt_nodes_only, boundary_node_id);  % nTrials × 3 logical
        dtr_BdNode(mask) = 1;

        %% constructing the big LMM/GLMM matrix
        group_col   = repmat(iExp, nTrials, 1);
        subj_col    = repmat(SubIdx, nTrials, 1);
        trialNo_col = log(trials_Col);
        dtNum_col   = dtNum - 1;
        cue_dtr_dis_col = cue_dtr_distance; % distance between cue and distractor
        tgt_dtr_dis_col = tgt_dtr_distance; % distance between target and distractor
        recency_col     = recency_nodes;
        walk_col        = rndHam_Col; % random or hamiltonian walk
        trans_col       = transStyle; % within or between transitions
        lureIn_col      = lureIn; 
        errFinal_col    = errorIdx_Col; % correct or incorrect based on the final click
        nodes_BdNode    = [cue_BdNode, tgt_BdNode, dtr_BdNode];
        RT_col          = respRT_Col;
        choice_col      = choiceId_ang;

        glmm_subj   = [glmm_subj; group_col, subj_col, trialNo_col, dtNum_col, ...
                                  cue_dtr_dis_col, tgt_dtr_dis_col, recency_col, ...
                                  walk_col, trans_col, lureIn_col, errFinal_col, ...
                                  nodes_BdNode, ...
                                  RT_col, choice_col];
    end
    %% save glmm_subj.mat
    save([bhvDataDir, ExpWord, '_glmm_subj.mat'], 'glmm_subj');

end


