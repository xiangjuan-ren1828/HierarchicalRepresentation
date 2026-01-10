% ClusterRep_v3_behvPaper_mouseTraj_biAssociation
% -----------------------------------------------
% modified based on ClusterRep_ChoiceModel_accuracy_association.m
% write by rxj @ 10/18/2022
%
% Can we found direct evidence to support association learning?
% --------Basic idea--------
% If participants experience A->B->C, then the accuracy should be higher if
% next target is B or A
% For a transition X->Y, the accuracy in the following 3 cases should be: (1) > (2) > (3)
% (1) Y occurs one trial before X;
% (2) Y occurs two trials before X;
% (3) other;
%
% --------Another calculation--------
% trial intervals of a transition occurs again: within vs. between
% transitions
% bidirectional or non-directional


clear
clc

%%
addpath('tight_subplot/');
addpath(genpath('HierarchicalCluster/'));

%% Columns index for the clusterResult.mat
rndTrial = 700;
HamTrial = 800;
nTrials  = rndTrial + HamTrial;

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

%%
circle_list = 0 : 1/60 : 1.5;
states      = 15;
nFit        = 100;
expList     = {'ImplicitExp', 'ExplicitExp', 'ImplicitRandExp'}; % {'ImplicitExp', 'ExplicitExp', 'ImplicitRandExp', 'ImplicitClusterExp', 'ImplicitRandClusterExp'};
folder      = '/Users/ren/Projects-NeuroCode/MyExperiment/HierarchicalCluster';

%%
subLen    = 24;
nBef      = 20;
nExp      = length(expList);
nTrans    = 7; %% within-trans, between-trans, random-walk+within-trans, random-walk+between-trans, hamiltonian-walk+within-trans, hamiltonian-walk+betwee-trans
transWord = {'within', 'between', 'rand-within', 'rand-between', 'ham-within', 'ham-between', 'all'};
acc_condition_group      = nan(subLen, nBef, nTrans, nExp); % nBef conditions
acc_nTgt_condition_group = nan(subLen, nBef, nTrans, nExp);
acc_rev_condition_group  = nan(subLen, nBef, nTrans, nExp); 
acc_rev_nTgt_condition_group = nan(subLen, nBef, nTrans, nExp);
tempInt_nonD_group = cell(subLen, nExp); %% only stimuli
tempInt_binD_group = cell(subLen, nExp); 
for iExp = 1 : nExp
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
    end
    subj_listBv = subj_list;
    subLen      = length(subj_list);
    
    %% loop across participants
    acc_condition_iSub          = nan(subLen, nBef, nTrans); % nBef conditions; 2: within and between transitions
    acc_nTgt_condition_iSub     = nan(subLen, nBef, nTrans); 
    acc_rev_condition_iSub      = nan(subLen, nBef, nTrans); % nBef conditions
    acc_rev_nTgt_condition_iSub = nan(subLen, nBef, nTrans); % nBef conditions
    for iSub = 1 : subLen
        disp([ExpWord, '-subj', num2str(iSub)]);
        %%
        subID  = subjLab{iSub};
        subjBv = subj_listBv{iSub};
        subjDir = [bhvDataDir, subjBv, '/'];
        load([subjDir, subjBv, 'clusterResult_Blocks.mat'], 'clusterResult');
        trials_Col   = (1 : 1 : nTrials)';
        respRT_Col   = clusterResult(:, respRT);
        errorIdx_Col = clusterResult(:, errorIdx);
        timeOIdx_Col = clusterResult(:, timeOIdx);
        rndHam_Col   = clusterResult(:, trialTp);
        %%% nodes identity
        from_nodes = clusterResult(:, objTgt);
        to_nodes   = clusterResult(:, objResp(1));
        dt_nodes   = clusterResult(:, objResp);
        %%% within vs. between cluster transition RTs across different distractor numbers
        dtNum = arrayfun(@(x) length(find(dt_nodes(x, :)~=0)), (1 : nTrials)');
        dtNum = dtNum - 1;
        %%% transition: within- or between-cluster transition
        transIn = zeros(nTrials, 3);
        transDt = zeros(nTrials, length(dtrNums)); % 3 kinds of distractors
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
            %%% for distractors
            for id = 1 : length(dtrNums)
                transDt(clsIdx, id) = arrayfun(@(x) ismember(x, clsDef(ic, :)), dt_nodes(clsIdx, id + 1));
            end
        end
        transStyle = (transIn(:, 1) == transIn(:, 2));
        
        %%% distance between cue image and distractor image
        cue_dtr_D = zeros(nTrials, length(dtrNums));
        for iT = 1 : nTrials
            from_i = from_nodes(iT);
            dt_i   = dt_nodes(iT, 2 : end);
            dt_i   = dt_i(dt_i ~= 0);
            cue_dtr_D(iT, 1 : length(dt_i)) = arrayfun(@(x) miniDmat(x, from_i), dt_i);
        end
        
        %%% mouse trajectory trials
        load([subjDir, subjBv,  'respYes_trials_blc.mat'], 'respYes_trials_blc');
        load([subjDir, subjBv,  'mouseTraj_trials_blc.mat'], 'mouseTraj_trials_blc');

        %%% extract the angles
        choiceId_ang  = nan(length(from_nodes), 1);
        dtrNoCnt_iSub = nan(length(from_nodes), 1);
        iCount = 1;
        for iBlock = 1 : nBlock
            blc_i = find((subNo_col == (iSub - 1)) & (blockNo_col == (iBlock - 1)));
            trlNo_i   = trialNo_col(blc_i);
            tgtAng_i  = tgtAng_col(blc_i);
            dtrAng1_i = dtrAng1_col(blc_i);
            dtrAng2_i = dtrAng2_col(blc_i);
            dtrAng3_i = dtrAng3_col(blc_i);
            dtrNoCnt_i= dtrNoCnt_col(blc_i);
            
            ttBlcIdx = (trialsInBlc(iBlock, 1) : trialsInBlc(iBlock, 2))';
            mouseTraj_trials = mouseTraj_trials_blc(ttBlcIdx);
            %load([subjDir, subjBv,  'mouseTraj_trials_in_Blc', num2str(iBlock), '.mat'], 'mouseTraj_trials');

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
                cueTp      = find(timePass <= 0.8);
                cueTp_stay = cueTp(end); %% the final stay point before cue onset
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
                iCount = iCount + 1;
            end
        end
        
        %% sliding across each transition pair
        %from_nodes = clusterResult(:, objTgt);
        %to_nodes   = clusterResult(:, objResp(1));
        acc_condition          = cell(nTrans, nBef); %% 2: within and between transitions; nBef conditions: (1) Y occurs one trial before X; (2) Y occurs two trials before X; (3) others;
        acc_nTgt_condition     = cell(nTrans, nBef); 
        acc_rev_condition      = cell(nTrans, nBef);
        acc_rev_nTgt_condition = cell(nTrans, nBef);
        for iT = 1 : nTrials
            from_i = from_nodes(iT);
            to_i   = to_nodes(iT);
            transStyle_i = transStyle(iT);
            rndHam_i     = rndHam_Col(iT);

            from_hty = from_nodes(1 : iT); %% history of the fromNodes
            to_hty   = to_nodes(1 : iT); 
            %% ---------------temporal association on choice---------------
            %% no matter whether to_i is transitioned to from_i
            tgtFind     = find(from_hty == to_i);
            tgtFind_sim = find(from_hty == from_i);
            tgtFind_sim(end) = [];
            nonFind = find(from_hty ~= to_i);
            nonFind(end) = [];
            %%% ----------For target----------
            if ~isempty(tgtFind)
                tgtFind_nearest = tgtFind(end);
                condIdx = length(from_hty) - tgtFind_nearest;
                if condIdx <= nBef 
                    if transStyle_i == 1     %% within-trans
                        acc_condition{1, condIdx} = [acc_condition{1, condIdx}; choiceId_ang(iT)];

                    elseif transStyle_i == 0 %% between-trans
                        acc_condition{2, condIdx} = [acc_condition{2, condIdx}; choiceId_ang(iT)];
                    end
                        
                    if transStyle_i == 1 && rndHam_i == 1 %% within-trans & random-walk
                        acc_condition{3, condIdx} = [acc_condition{3, condIdx}; choiceId_ang(iT)];

                    elseif transStyle_i == 0 && rndHam_i == 1 %% between-trans & random-walk
                        acc_condition{4, condIdx} = [acc_condition{4, condIdx}; choiceId_ang(iT)];

                    elseif transStyle_i == 1 && rndHam_i == 2 %% within-trans & hamiltonian-walk
                        acc_condition{5, condIdx} = [acc_condition{5, condIdx}; choiceId_ang(iT)];

                    elseif transStyle_i == 0 && rndHam_i == 2 %% between-trans & hamiltonian-walk
                        acc_condition{6, condIdx} = [acc_condition{6, condIdx}; choiceId_ang(iT)];
                    end
                    acc_condition{7, condIdx} = [acc_condition{7, condIdx}; choiceId_ang(iT)];
                end
            end
            %%% ----------For non-target----------
            %%% Is it possible the condIdx larger than 1???????
            if ~isempty(nonFind) 
                nonFind_nearest = nonFind(end);
                if (isempty(tgtFind_sim) || (tgtFind_sim(end) < nonFind_nearest)) && ...
                   (isempty(tgtFind) || (tgtFind(end) < nonFind_nearest))
                    condIdx = length(from_hty) - nonFind_nearest;
                    if condIdx <= nBef
                        if transStyle_i == 1     %% within-trans
                            acc_nTgt_condition{1, condIdx} = [acc_nTgt_condition{1, condIdx}; choiceId_ang(iT)];

                        elseif transStyle_i == 0 %% between-trans
                            acc_nTgt_condition{2, condIdx} = [acc_nTgt_condition{2, condIdx}; choiceId_ang(iT)];
                        end

                        if transStyle_i == 1 && rndHam_i == 1     %% within-trans & random-walk
                            acc_nTgt_condition{3, condIdx} = [acc_nTgt_condition{3, condIdx}; choiceId_ang(iT)];

                        elseif transStyle_i == 0 && rndHam_i == 1 %% between-trans & random-walk
                            acc_nTgt_condition{4, condIdx} = [acc_nTgt_condition{4, condIdx}; choiceId_ang(iT)];

                        elseif transStyle_i == 1 && rndHam_i == 2 %% within-trans & hamiltonian-walk
                            acc_nTgt_condition{5, condIdx} = [acc_nTgt_condition{5, condIdx}; choiceId_ang(iT)];

                        elseif transStyle_i == 0 && rndHam_i == 2 %% between-trans & hamiltonian-walk
                            acc_nTgt_condition{6, condIdx} = [acc_nTgt_condition{6, condIdx}; choiceId_ang(iT)];
                        end
                        acc_nTgt_condition{7, condIdx} = [acc_nTgt_condition{7, condIdx}; choiceId_ang(iT)];
                    end
                end
            end

            %% to_i occurs before and associated with from_i
            %toFind_ass = find((from_hty == from_i & to_hty == to_i) | (from_hty == to_i & to_hty == from_i));
            %toFind_ass(end) = []; % remove the self finding
            tgtFind_ass = find((from_hty == to_i & to_hty == from_i));
            tgtFind_sim = find((from_hty == from_i & to_hty == to_i));
            tgtFind_sim(end) = [];
            nonFind_ass = find(~(from_hty == to_i & to_hty == from_i)); %find((from_hty ~= to_i & to_hty ~= from_i));
            nonFind_ass(end) = [];
            %%% ----------For target----------
            if ~isempty(tgtFind_ass)
                tgtFind_ass_nearest = tgtFind_ass(end);
                if isempty(tgtFind_sim) || (tgtFind_sim(end) < tgtFind_ass_nearest)
                    condIdx = length(from_hty) - tgtFind_ass_nearest;
                    if condIdx <= nBef
                        if transStyle_i == 1    %% within-trans 
                            acc_rev_condition{1, condIdx} = [acc_rev_condition{1, condIdx}; choiceId_ang(iT)];

                        elseif transStyle_i == 0 %% between-trans 
                            acc_rev_condition{2, condIdx} = [acc_rev_condition{2, condIdx}; choiceId_ang(iT)];
                        end

                        if transStyle_i == 1 && rndHam_i == 1     %% within-trans & random-walk
                            acc_rev_condition{3, condIdx} = [acc_rev_condition{3, condIdx}; choiceId_ang(iT)];

                        elseif transStyle_i == 0 && rndHam_i == 1 %% between-trans & random-walk
                            acc_rev_condition{4, condIdx} = [acc_rev_condition{4, condIdx}; choiceId_ang(iT)];

                        elseif transStyle_i == 1 && rndHam_i == 2 %% within-trans & hamiltonian-walk
                            acc_rev_condition{5, condIdx} = [acc_rev_condition{5, condIdx}; choiceId_ang(iT)];

                        elseif transStyle_i == 0 && rndHam_i == 2 %% between-trans & hamiltonian-walk
                            acc_rev_condition{6, condIdx} = [acc_rev_condition{6, condIdx}; choiceId_ang(iT)];
                        end
                        acc_rev_condition{7, condIdx} = [acc_rev_condition{7, condIdx}; choiceId_ang(iT)]; % not differentiate 4 transition conditions
                    end

                end
            end
            %%% ----------For non-target----------
            if ~isempty(nonFind_ass)
                nonFind_ass_nearest = nonFind_ass(end);
                if (isempty(tgtFind_sim) || (tgtFind_sim(end) < nonFind_ass_nearest)) && ...
                   (isempty(tgtFind_ass) || (tgtFind_ass(end) < nonFind_ass_nearest))
                    condIdx = length(from_hty) - nonFind_ass_nearest;
                    if condIdx <= nBef 
                        if transStyle_i == 1     %% within-trans & random-walk
                            acc_rev_nTgt_condition{1, condIdx} = [acc_rev_nTgt_condition{1, condIdx}; choiceId_ang(iT)];

                        elseif transStyle_i == 0 %% between-trans & random-walk
                            acc_rev_nTgt_condition{2, condIdx} = [acc_rev_nTgt_condition{2, condIdx}; choiceId_ang(iT)];
                        end

                        if transStyle_i == 1 && rndHam_i == 1     %% within-trans & random-walk
                            acc_rev_nTgt_condition{3, condIdx} = [acc_rev_nTgt_condition{3, condIdx}; choiceId_ang(iT)];

                        elseif transStyle_i == 0 && rndHam_i == 1 %% between-trans & random-walk
                            acc_rev_nTgt_condition{4, condIdx} = [acc_rev_nTgt_condition{4, condIdx}; choiceId_ang(iT)];

                        elseif transStyle_i == 1 && rndHam_i == 2 %% within-trans & hamiltonian-walk
                            acc_rev_nTgt_condition{5, condIdx} = [acc_rev_nTgt_condition{5, condIdx}; choiceId_ang(iT)];

                        elseif transStyle_i == 0 && rndHam_i == 2 %% between-trans & hamiltonian-walk
                            acc_rev_nTgt_condition{6, condIdx} = [acc_rev_nTgt_condition{6, condIdx}; choiceId_ang(iT)];
                        end
                        acc_rev_nTgt_condition{7, condIdx} = [acc_rev_nTgt_condition{7, condIdx}; choiceId_ang(iT)]; % not differentiate 4 transition conditions
                    end
                end
            end

            %% ---------------temporal intervals of last occurrence---------------
            % When was the last time the transition pair appeared
            % bidirection vs. non-direction
            tgtFind_nonD = find((from_hty == from_i & to_hty == to_i) | (from_hty == to_i & to_hty == from_i)); % nondirectional
            tgtFind_binD = find(from_hty == from_i & to_hty == to_i); % directional
            tgtFind_nonD(end) = [];
            tgtFind_binD(end) = [];
            if ~isempty(tgtFind_nonD)
                tgtFind_nonD_nearest = tgtFind_nonD(end);
                condIdx = length(from_hty) - tgtFind_nonD_nearest;

                if transStyle_i == 1 && rndHam_i == 1     %% within-trans & random-walk
                    tempInt_i = [condIdx, 1, iT];

                elseif transStyle_i == 0 && rndHam_i == 1 %% between-trans & random-walk
                    tempInt_i = [condIdx, 2, iT];

                elseif transStyle_i == 1 && rndHam_i == 2 %% within-trans & hamiltonian-walk
                    tempInt_i = [condIdx, 3, iT];

                elseif transStyle_i == 0 && rndHam_i == 2 %% between-trans & hamiltonian-walk
                    tempInt_i = [condIdx, 4, iT];
                end
                tempInt_nonD_group{iSub, iExp} = [tempInt_nonD_group{iSub, iExp}; tempInt_i];
            end

            if ~isempty(tgtFind_binD)
                tgtFind_binD_nearest = tgtFind_binD(end);
                condIdx = length(from_hty) - tgtFind_binD_nearest;
                if transStyle_i == 1 && rndHam_i == 1     %% within-trans & random-walk
                    tempInt_i = [condIdx, 1, iT];

                elseif transStyle_i == 0 && rndHam_i == 1 %% between-trans & random-walk
                    tempInt_i = [condIdx, 2, iT];

                elseif transStyle_i == 1 && rndHam_i == 2 %% within-trans & hamiltonian-walk
                    tempInt_i = [condIdx, 3, iT];

                elseif transStyle_i == 0 && rndHam_i == 2 %% between-trans & hamiltonian-walk
                    tempInt_i = [condIdx, 4, iT];
                end
                tempInt_binD_group{iSub, iExp} = [tempInt_binD_group{iSub, iExp}; tempInt_i];
            end

        end
        for iCon = 1 : nBef
            for iTrans = 1 : nTrans
                acc_iCon = acc_condition{iTrans, iCon};
                acc_iCon(isnan(acc_iCon)) = [];
                acc_condition_iSub(iSub, iCon, iTrans) = length(find(acc_iCon == 1)) / length(acc_iCon);

                acc_nTgt_iCon = acc_nTgt_condition{iTrans, iCon};
                acc_nTgt_iCon(isnan(acc_nTgt_iCon)) = [];
                acc_nTgt_condition_iSub(iSub, iCon, iTrans) = length(find(acc_nTgt_iCon == 1)) / length(acc_nTgt_iCon); 

                acc_rev_iCon = acc_rev_condition{iTrans, iCon};
                acc_rev_iCon(isnan(acc_rev_iCon)) = [];
                acc_rev_condition_iSub(iSub, iCon, iTrans) = length(find(acc_rev_iCon == 1)) / length(acc_rev_iCon);

                acc_rev_nTgt_iCon = acc_rev_nTgt_condition{iTrans, iCon};
                acc_rev_nTgt_iCon(isnan(acc_rev_nTgt_iCon)) = [];
                acc_rev_nTgt_condition_iSub(iSub, iCon, iTrans) = length(find(acc_rev_nTgt_iCon == 1)) / length(acc_rev_nTgt_iCon);

            end
        end
    end
    acc_condition_group(:, :, :, iExp)      = acc_condition_iSub;
    acc_nTgt_condition_group(:, :, :, iExp) = acc_nTgt_condition_iSub;
    acc_rev_condition_group(:, :, :, iExp)  = acc_rev_condition_iSub;
    acc_rev_nTgt_condition_group(:, :, :, iExp) = acc_rev_nTgt_condition_iSub;

end
%% save 'acc_condition_group' and 'acc_nTgt_condition_group'
% save([folder, '/acc_condition_group-TrueData.mat'], 'acc_condition_group');
% save([folder, '/acc_nTgt_condition_group-TrueData.mat'], 'acc_nTgt_condition_group');

%% color settings
colorSets = [0.98, 0.72, 0.69; ...
            0.97, 0.85, 0.67; ...
            0.33, 0.73, 0.83; ...
            0.72, 0.80, 0.88; ...
            0.54, 0.67, 0.20; ...
            0.82, 0.92, 0.78; ...
            0.78, 0.50, 0.75; ...
            0.86, 0.80, 0.89; ...
            0.75, 0.56, 0; ...
            0.40, 0.40, 0.40];
colorComp = colorSets([1, 2, 3, 4, 7], :);

%% plotting
% acc_condition_group      = nan(subLen, nBef, nTrans, nExp); % nBef conditions
% acc_nTgt_condition_group = nan(subLen, nBef, nTrans, nExp);
% acc_rev_condition_group  = nan(subLen, nBef, nTrans, nExp); 
% acc_rev_nTgt_condition_group = nan(subLen, nBef, nTrans, nExp);
accFlg = 2;
if accFlg == 1
    acc_plot = acc_condition_group;
    non_plot = acc_nTgt_condition_group;
elseif accFlg == 2
    acc_plot = acc_rev_condition_group;
%     acc_plot(:, 2, :, 1) = nan;
    non_plot = acc_rev_nTgt_condition_group;
end

%% SI figure: plotting target and non-target in each condition
for iExp = 3% : nExp
    acc_plot_iExp = acc_plot(:, :, :, iExp);
    non_plot_iExp = non_plot(:, :, :, iExp);
    [acc_avg, acc_sem] = Mean_and_Se(acc_plot_iExp, 1);
    acc_avg = squeeze(acc_avg); % nBef * nTrans
    acc_sem = squeeze(acc_sem);
    [non_avg, non_sem] = Mean_and_Se(non_plot_iExp, 1);
    non_avg = squeeze(non_avg); 
    non_sem = squeeze(non_sem);
    for iTrans = 1 : nTrans % within-trans and between-trans
        figure('Position', [100 100 380 240]), clf;
        % target
        errorbar(1 : 1 : nBef, acc_avg(:, iTrans), acc_sem(:, iTrans), 'Color', colorComp(iTrans, :), 'LineStyle', '-', 'LineWidth', 2); hold on;
        plot(1 : 1 : nBef, acc_avg(:, iTrans), 'Marker', '.', 'MarkerSize', 20, 'Color', colorComp(iTrans, :), 'LineStyle', 'none'); hold on;
        % non-target
        errorbar(1 : 1 : nBef, non_avg(:, iTrans), non_sem(:, iTrans), 'Color', [0.4, 0.4, 0.4], 'LineStyle', '-', 'LineWidth', 2); hold on;
        plot(1 : 1 : nBef, non_avg(:, iTrans), 'Marker', '.', 'MarkerSize', 20, 'Color', [0.4, 0.4, 0.4], 'LineStyle', 'none'); hold on;
        %plot(xlim, [0.5, 0.5], 'k--', 'LineWidth', 1); hold on;
        plot(xlim, [(1/2+1/3+1/4)/3, (1/2+1/3+1/4)/3], 'k--', 'LineWidth', 1); hold on;
        ylim([0, 1]);
        set(gca, 'LineWidth', 2);
        set(gca, 'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Arial');
        box off;
        %%% statistical tests
        disp(['======', expList{iExp}, ', ', transWord{iTrans}, '======'])
        [h, p, ci, stats] = ttest(acc_plot_iExp(:, 1, iTrans), non_plot_iExp(:, 1, iTrans))
    end
    %ylim([0, 1]);
    %ylimit = ylim;
end

%% BehaviorPaper Figure xx: within vs. between transitions, bi-directional associations, if X-Y happens in the current trial, whether Y-X/non Y-X in the laster trial has a influence
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

% acc_trans_walk = zeros(subLen, 2, 2, length(expList)); % the 1st and 2nd '2' denote 'within vs. between cluster transition' and 'random and hamiltonian' walk
barPos   = [1, 1.5; 1.7, 2.2];
chance_i = mean([1/2, 1/3, 1/4]);

iPrev = 1;
acc_plot_nExp = squeeze(acc_plot(:, iPrev, 1 : 2, :)); % nan(subLen, nTrans, nExp); nTrans: within vs. between transition
non_plot_nExp = squeeze(non_plot(:, iPrev, 1 : 2, :));
for iExp = 1 : nExp
    disp(['---------- ', expList{iExp}, ' ----------']);
    figure('Position', [100 100 260 120]), clf;

    for i_tp = 1 : 2 % within vs. between-transition
        if i_tp == 1
            transWord = 'Within-trans';
        elseif i_tp == 2
            transWord = 'Between-trans';
        end
        barPos_i = barPos(i_tp, :);

        ass_plot = nan(size(acc_plot_nExp, 1), 2); % with the reversed directional transition happening vs. not
        ass_plot(:, 1) = acc_plot_nExp(:, i_tp, iExp);
        ass_plot(:, 2) = non_plot_nExp(:, i_tp, iExp);
        [assAvg, assSem] = Mean_and_Se(ass_plot, 1);

        %%% line plot
        plot(barPos_i, ass_plot, 'Color', [0.6, 0.6, 0.6], 'LineStyle', '-', 'LineWidth', indvLineW); hold on;
        plot(barPos_i, assAvg, 'Color', [0, 0, 0], 'LineStyle', '-', 'LineWidth', errLineWid); hold on;
        for i_pre = 1 : 2 % with the reversed directional transition happening vs. not
            if i_pre == 1
                preWord  = 'Y-X before';
                colorTmp = [0, 0, 0];
            elseif i_pre == 2
                preWord  = 'Y-X none';
                colorTmp = [1, 1, 1];
            end
            errorbar(barPos_i(i_pre), assAvg(i_pre), assSem(i_pre), 'Color', 'k', 'LineStyle', 'none', 'LineWidth', errLineWid); hold on;
            plot(barPos_i(i_pre), assAvg(i_pre), 'Marker', 'o', 'MarkerSize', markSize, 'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', colorTmp, 'LineStyle', '-'); hold on;
        end
        % ------ Statistical tests ------
        disp(['======== ', transWord, '-Y-X before vs. Y-X none ========']);
        [h, p, ci, stats] = ttest(ass_plot(:, 1), ass_plot(:, 2))
    end
    xlim([0.6, 2.6]);
    ylim([0, 1]);
    plot(xlim, [chance_i, chance_i], 'k--', 'LineWidth', 0.8); hold on;
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
    save_name = ['Exp', num2str(iExp), '-', expList{iExp}, '-biAss-trans.png'];
    exportgraphics(ax, save_name, 'Resolution', 600);

end

%% SI figure: within vs. between transitions in Random and Hamiltonian Walks, bi-directional associations, if X-Y happens in the current trial, whether Y-X/non Y-X in the laster trial has a influence


