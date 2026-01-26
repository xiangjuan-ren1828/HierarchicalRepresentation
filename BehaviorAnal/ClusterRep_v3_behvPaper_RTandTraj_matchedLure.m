% ClusterRep_v3_behvPaper_RTandTraj_matchedLure.m
% modified by XR @ Dec 30 2025 based on the followings
% 
% 1) RT or choiceAccuracy with lure distractors –– matched the distance
% between the cue to the lure stimulus or to the non-lure stimulus
% 2) For the error trials with both lure and another 2-step control
% stimulus as the distractor, whether participants have a larger
% probability to choose the lure over the control
% =========================================================================
% ClusterRep_v3_RTAnal_lureDtr
% write by rxj @ 10/25/2021
%

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
rndInter = 85;
HamInter = 15;
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


%% For the boundary node, if another boundary within the same cluster exists, how about the RTs?
% added by rxj @ 08/24/2021
% !!!!!! match the distance between the lure stimulus and non-lure stimulus !!!!!!
%% within-cluster transtions
bndNode_Num = 6;
transIn_clsInDtr = cell(1, bndNode_Num);
transIn_clsInDtr{1} = [1, 2, 5; 1, 3, 5; 1, 4, 5]; % fromNode-toNode-dtrNode
transIn_clsInDtr{2} = [5, 2, 1; 5, 3, 1; 5, 4, 1];
transIn_clsInDtr{3} = [6, 7, 10; 6, 8, 10; 6, 9, 10];
transIn_clsInDtr{4} = [10, 7, 6; 10, 8, 6; 10, 9, 6];
transIn_clsInDtr{5} = [11, 12, 15; 11, 13, 15; 11, 14, 15];
transIn_clsInDtr{6} = [15, 12, 11; 15, 13, 11; 15, 14, 11];

transIn_clsOutDtr = cell(1, bndNode_Num);
transIn_clsOutDtr{1} = [1, 2, 12; ...
                        1, 2, 13; ...
                        1, 2, 14; ...
                        1, 3, 12; ...
                        1, 3, 13; ...
                        1, 3, 14; ...
                        1, 4, 12; ...
                        1, 4, 13; ...
                        1, 4, 14];
transIn_clsOutDtr{2} = [5, 2, 7; ...
                        5, 2, 8; ...
                        5, 2, 9; ...
                        5, 3, 7; ...
                        5, 3, 8; ...
                        5, 3, 9; ...
                        5, 4, 7; ...
                        5, 4, 8; ...
                        5, 4, 9];
transIn_clsOutDtr{3} = [6, 7, 2; ...
                        6, 7, 3; ...
                        6, 7, 4; ...
                        6, 8, 2; ...
                        6, 8, 3; ...
                        6, 8, 4; ...
                        6, 9, 2; ...
                        6, 9, 3; ...
                        6, 9, 4];
transIn_clsOutDtr{4} = [10, 7, 12; ...
                        10, 7, 13; ...
                        10, 7, 14; ...
                        10, 8, 12; ...
                        10, 8, 13; ...
                        10, 8, 14; ...
                        10, 9, 12; ...
                        10, 9, 13; ...
                        10, 9, 14];
transIn_clsOutDtr{5} = [11, 12, 7; ...
                        11, 12, 8; ...
                        11, 12, 9; ...
                        11, 13, 7; ...
                        11, 13, 8; ...
                        11, 13, 9; ...
                        11, 14, 7; ...
                        11, 14, 8; ...
                        11, 14, 9];
transIn_clsOutDtr{6} = [15, 12, 2; ...
                        15, 12, 3; ...
                        15, 12, 4; ...
                        15, 13, 2; ...
                        15, 13, 3; ...
                        15, 13, 4; ....
                        15, 14, 2; ...
                        15, 14, 3; ...
                        15, 14, 4];

%% between-cluster transtions
transOut_clsInDtr = cell(1, bndNode_Num);
transOut_clsInDtr{1} = [1, 15, 5];
transOut_clsInDtr{2} = [5, 6, 1];
transOut_clsInDtr{3} = [6, 5, 10];
transOut_clsInDtr{4} = [10, 11, 6];
transOut_clsInDtr{5} = [11, 10, 15];
transOut_clsInDtr{6} = [15, 1, 11];

transOut_clsOutDtr = cell(1, bndNode_Num);
transOut_clsOutDtr{1} = [1, 15, 12; 1, 15, 13; 1, 15, 14];
transOut_clsOutDtr{2} = [5, 6, 7; 5, 6, 8; 5, 6, 9];
transOut_clsOutDtr{3} = [6, 5, 2; 6, 5, 3; 6, 5, 4];
transOut_clsOutDtr{4} = [10, 11, 12; 10, 11, 13; 10, 11, 14];
transOut_clsOutDtr{5} = [11, 10, 7; 11, 10, 8; 11, 10, 9];
transOut_clsOutDtr{6} = [15, 1, 2; 15, 1, 3; 15, 1, 4];

%% indicating the distractor
dtrCond = input('DistractorNo condition, 1-one, 2-two, 3-three, 4-merged: '); % only trials with 1, 2, 3 distractors or 4 (merge all kinds of distractors)

%% 3 Experiments:
% I didn't regress out the confounding variables from RTs.
circle_list = 0 : 1/60 : 1.5;
states      = 15;
nFit        = 100;
expList     = {'ImplicitExp', 'ExplicitExp', 'ImplicitRandExp'}; 
folder      = '/Users/ren/Projects-NeuroCode/MyExperiment/HierarchicalCluster';

subLen = 24;
% ---------- RTs ----------
RT_transIn_subj        = nan(bndNode_Num, 2, subLen, length(expList)); % 2: with vs. without lure stimulus
RT_transOut_subj       = nan(bndNode_Num, 2, subLen, length(expList));
transIn_Len            = nan(bndNode_Num, 2, subLen, length(expList));
transOut_Len           = nan(bndNode_Num, 2, subLen, length(expList));
% ---------- Choice accuracy: one-hot measure ----------
acc_transIn_subj       = nan(bndNode_Num, 2, subLen, length(expList));
acc_transOut_subj      = nan(bndNode_Num, 2, subLen, length(expList));
% ---------- Choice accuracy: trajectory ----------
acc_transIn_traj_subj  = nan(bndNode_Num, length(circle_list), 2, subLen, length(expList));
acc_transOut_traj_subj = nan(bndNode_Num, length(circle_list), 2, subLen, length(expList));
% ---------- Quantify the distance between cue and distractors, distance between target and distractors ----------
cue_dtr_distance_subj  = cell(subLen, 4, length(expList)); % 4: with vs. without lure stimulus for boundary-to-within and boundary-to-boundary transitions
tgt_dtr_distance_subj  = cell(subLen, 4, length(expList));
% ---------- Relative angles between the mouse and the target/lure-distractor/control ----------
ang_transIn_subj       = nan(bndNode_Num, 2, subLen, length(expList));
ang_transOut_subj      = nan(bndNode_Num, 2, subLen, length(expList));

% ---------- In the correct response trials with boundary node transition
% and lure + 2-step control stimulus, if participants have a larger
% probability to select the lure >> control >> others ----------
lureEffect_errTrial_transIn_subj  = nan(bndNode_Num, 3, subLen, length(explist));
lureEffect_errTrial_transOut_subj = nan(bndNode_Num, 3, subLen, length(explist));

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

    for SubIdx = 1 : subLen
        disp([ExpWord, '-subj', num2str(SubIdx)]);
        %%
        subjBv      = subj_listBv{SubIdx};
        subTit      = subjLab{SubIdx};
        subjDir     = [bhvDataDir, subjBv, '/'];
        load([subjDir, subjBv, 'clusterResult_Blocks.mat'], 'clusterResult');

        respRT_Col   = clusterResult(:, respRT);
        errorIdx_Col = clusterResult(:, errorIdx);
        timeOIdx_Col = clusterResult(:, timeOIdx);
        rndHam_Col   = clusterResult(:, trialTp);
        fast_Col     = clusterResult(:, fastIdx);
        %%% within vs. between cluster transition RTs
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
        transStyle = (transIn(:, 1) == transIn(:, 2)); % 1-within transition; 0-between transition
        nodesLabel = zeros(nTrials, 1);
        nodesLabel(from_nodes == 1 | from_nodes == 5 | from_nodes == 6 | from_nodes == 10 | from_nodes == 11 | from_nodes == 15) = 1;

        %%% within vs. between cluster transition RTs across different distractor numbers
        dtNum = arrayfun(@(x) length(find(dt_nodes(x, :)~=0)), (1 : nTrials)');
        dtNum = dtNum - 1;

        %%% mouse trajectory trials
        load([subjDir, subjBv,  'respYes_trials_blc.mat'], 'respYes_trials_blc');
        load([subjDir, subjBv,  'mouseTraj_trials_blc.mat'], 'mouseTraj_trials_blc');

        %%% !!!!!!!! remove the time-out trials !!!!!!!! %%%
        stim = clusterResult(:, objTgt);
        resp = clusterResult(:, objResp);
        %%% extract the angles
        choiceId_ang  = nan(length(stim), 1);
        choiceId      = nan(length(stim), length(circle_list));
        dtrNoCnt_iSub = nan(length(stim), 1);
        %%% angle for the target and all distractors
        angleMag_oneH = nan(length(stim), 4);
        angleMag_time = nan(length(stim), length(circle_list), 4);

        %%% Determine the identity of the chosen option
        chosen_opt_trials = nan(length(stim), 1); % the column index of which object was chosen in dt_nodes(iCount, :)
        iCount   = 1;
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
            %load([subjDir, subjBv,  'mouseTraj_trials_in_Blc', num2str(iBlock), '.mat'], 'mouseTraj_trials');

            for jTrl = 1 : length(mouseTraj_trials)
                trl_j     = find(trlNo_i == (jTrl - 1));
                tgtAng_j  = tgtAng_i(trl_j);
                dtrAng1_j = dtrAng1_i(trl_j);
                dtrAng2_j = dtrAng2_i(trl_j);
                dtrAng3_j = dtrAng3_i(trl_j);
                objAngs   = [tgtAng_j, dtrAng1_j, dtrAng2_j, dtrAng3_j];
                dtrNoCnt_j= dtrNoCnt_i(trl_j);
                dtrNoCnt_iSub(iCount)    = unique(dtrNoCnt_j);
   
                timePass  = mouseTraj_trials{jTrl}(:, 4);
                %%% click time point (response before cue) or time point
                %%% before cue onset (response after cue)
                cueTp             = find(timePass <= 0.8);
                cueTp_stay        = cueTp(end); %% the final stay point before cue onset
                objAngs_stay      = objAngs(cueTp_stay, :);
                objAngs_stay_left = objAngs_stay(1 : (unique(dtrNoCnt_j) + 1));
                if ~any(isnan(objAngs_stay_left))
                    [~, minId] = min(abs(objAngs_stay_left));
                    if minId == 1
                        choiceId_ang(iCount) = 1; % 1: choose the target

                    elseif minId ~= 1
                        choiceId_ang(iCount) = 0; % 0: choose the non-target
                    end
                    angleMag_oneH(iCount, 1 : length(objAngs_stay_left)) = abs(objAngs_stay_left);
                    % ------ which option participant chose ------
                    chosen_opt_trials(iCount) = minId;
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
                            angleMag_time(iCount, iTp, 1 : length(objAngs_stay_left)) = abs(objAngs_stay_left);
                        end
                    end
                end
                iCount = iCount + 1;
            end
        end

        % ------------ distractor conditions ------------
        if dtrCond == 1 || dtrCond == 2 || dtrCond == 3
            errorIdx_Col = errorIdx_Col(dtNum == dtrCond, :);
            fast_Col     = fast_Col(dtNum == dtrCond, :);
            rndHam_Col   = rndHam_Col(dtNum == dtrCond, :);
            from_nodes   = from_nodes(dtNum == dtrCond, :);
            to_nodes     = to_nodes(dtNum == dtrCond, :);
            dt_nodes     = dt_nodes(dtNum == dtrCond, 2 : (dtrCond + 1));
            respRT_Col   = respRT_Col(dtNum == dtrCond, :);
            choiceId_ang = choiceId_ang(dtNum == dtrCond, :);
            choiceId     = choiceId(dtNum == dtrCond, :);
        else
            dt_nodes     = dt_nodes(:, 2 : end);
        end

        %% quantify the influence of lure stimulus on RTs and choice accuracy
        for i = 1 : bndNode_Num
            %% ---------transition: from boundary node to within node---------
            transIn_InDtr  = transIn_clsInDtr{i}; % transIn_clsInDtr{1} = [1, 2, 5; 1, 3, 5; 1, 4, 5];
            transIn_OutDtr = transIn_clsOutDtr{i};
            exclude_InDtr  = unique(transIn_OutDtr(:, 3)); % should be 3 nodes
            exclude_OutDtr = unique(transIn_InDtr(:, 3));  % should be 1 node
            % **************** 1. boundary-to-within transition: with lure ****************
            RTs_InDtr = [];
            acc_InDtr = [];
            acc_InDtr_traj = cell(1, length(circle_list));
            ang_InDtr      = [];
            ang_InDtr_traj = cell(1, length(circle_list));
            for j = 1 : size(transIn_InDtr, 1)
                pair_j = transIn_InDtr(j, :);
                if dtrCond == 1
                    dt_Yes = sum((dt_nodes == pair_j(3)), 2);
                else
                    no_control = arrayfun(@(x) ~any(ismember(dt_nodes(x, :), exclude_InDtr)), (1 : length(dt_nodes))');
                    dt_Yes = sum((dt_nodes == pair_j(3)), 2) & no_control;
                end
                % ------ RTs ------
                RTs_InDtr = [RTs_InDtr; respRT_Col(find(from_nodes == pair_j(1) & to_nodes == pair_j(2) & dt_Yes == 1 & errorIdx_Col ~= 1 & fast_Col ~= 1))];
                % ------ choice accuracy ------
                acc_InDtr = [acc_InDtr; length(find(choiceId_ang == 1 & from_nodes == pair_j(1) & to_nodes == pair_j(2) & dt_Yes == 1)), length(find(~isnan(choiceId_ang) & from_nodes == pair_j(1) & to_nodes == pair_j(2) & dt_Yes == 1))];
                % ------ choice accuracy trajectory ------
                for iTp = 1 : length(circle_list)
                    choiceId_iTp = choiceId(:, iTp);
                    acc_InDtr_traj{iTp} = [acc_InDtr_traj{iTp}; length(find(choiceId_iTp == 1 & from_nodes == pair_j(1) & to_nodes == pair_j(2) & dt_Yes == 1)), length(find(~isnan(choiceId_iTp) & from_nodes == pair_j(1) & to_nodes == pair_j(2) & dt_Yes == 1))];
                end
                % ------ quantify the distance between cue/target and the non-lure distractor ------
                cue_j = pair_j(1);
                tgt_j = pair_j(2);
                dtr_j = dt_nodes(from_nodes == pair_j(1) & to_nodes == pair_j(2) & dt_Yes == 1, :);
                for dd = 1 : size(dtr_j, 1)
                    dtr_j_dd = dtr_j(dd, :);
                    dtr_j_dd = dtr_j_dd(dtr_j_dd ~= 0 & dtr_j_dd ~= pair_j(3));
                    minD_cue_i = miniDmat(dtr_j_dd', cue_j);
                    minD_tgt_i = miniDmat(dtr_j_dd', tgt_j);

                    cue_dtr_distance_subj{SubIdx, 1, iExp} = [cue_dtr_distance_subj{SubIdx, 1, iExp}; minD_cue_i];
                    tgt_dtr_distance_subj{SubIdx, 1, iExp} = [tgt_dtr_distance_subj{SubIdx, 1, iExp}; minD_tgt_i];
                end
                % ------ quantify the relative radian between ang(mouse-target) and ang(mouse-lure/control) ------
                ang_tgt      = angleMag_oneH(find(from_nodes == pair_j(1) & to_nodes == pair_j(2) & dt_Yes == 1), 1);
                ang_lure_tmp = angleMag_oneH(find(from_nodes == pair_j(1) & to_nodes == pair_j(2) & dt_Yes == 1), 2 : end);
                dt_nodes_tmp = dt_nodes(find(from_nodes == pair_j(1) & to_nodes == pair_j(2) & dt_Yes == 1), :);
                ang_lure_tmp = ang_lure_tmp(~isnan(ang_tgt), :);
                dt_nodes_tmp = dt_nodes_tmp(~isnan(ang_tgt), :);
                ang_tgt      = ang_tgt(~isnan(ang_tgt));
                ang_lure     = nan(size(ang_lure_tmp, 1), 1);
                for ll = 1 : size(ang_lure_tmp, 1)
                    ang_lure_tmp_ll = ang_lure_tmp(ll, :);
                    ang_lure(ll) = ang_lure_tmp_ll(find(dt_nodes_tmp(ll, :) == pair_j(3)));
                end
                ang_InDtr = [ang_InDtr; ang_tgt, ang_lure];

                % ---------- Competition between the Lure stimulus and the 2-step control stimulus ---------- 
                control_stim_mat  = arrayfun(@(x) ismember(dt_nodes(x, :), exclude_InDtr), (1 : length(dt_nodes))'); % if control stimulus exists
                lure_stim_mat     = (dt_nodes == pair_j(3)); % if lure stimulus exists
                % ------ Get the column index with 1 for each row ------
                % For 'control_stim_mat', there might be two '1' for specific trials (rows)
                [row_idx, col_idx] = find(control_stim_mat);
                % Put results into a cell array: one cell per row
                control_stim_idx   = accumarray(row_idx, col_idx, [], @(x){x});
                [~, lure_stim_idx] = max(lure_stim_mat, [], 2);



                control_stim_yes = arrayfun(@(x) any(ismember(dt_nodes(x, :), exclude_InDtr)), (1 : length(dt_nodes))');
                lure_stim_yes    = sum((dt_nodes == pair_j(3)), 2);
                error_choice_yes = choiceId_ang == 0;
                comp_trials_yes  = (control_stim_yes) & (lure_stim_yes) & (error_choice_yes);
                if sum(comp_trials_yes) ~= 0
                    % proportation that participant chose the lure, the
                    % control and the other
                    comp_idx = find(comp_trials_yes == 1);
                    for i_cm = 1 : length(comp_idx)
                        chosen_opt_ii = chosen_opt_trials(comp_idx(i_cm));
                        if ~isnan(chosen_opt_ii)
                            if chosen_opt_ii == lure_stim_idx(comp_idx(i_cm)) % choosing lure stimulus

                            elseif ismember(chosen_opt_ii, control_stim_idx{comp_idx(i_cm)}) % choosing the control
                                
                            else % choosing the other

                            end
                        end
                    end

                end
                


            end
            if ~isempty(ang_InDtr)
                ang_InDtr_diff = mod(ang_InDtr(:, 2) - ang_InDtr(:, 1) + pi, 2*pi) - pi;
                ang_transIn_subj(i, 1, SubIdx, iExp) = atan2(mean(sin(ang_InDtr_diff)), mean(cos(ang_InDtr_diff)));
            end

            RT_transIn_subj(i, 1, SubIdx, iExp)  = nanmean(RTs_InDtr);
            acc_transIn_subj(i, 1, SubIdx, iExp) = sum(acc_InDtr(:, 1)) / sum(acc_InDtr(:, 2));

            % **************** 2. boundary-to-within transition: without lure ****************
            RTs_OutDtr = [];
            acc_OutDtr = [];
            acc_OutDtr_traj = cell(1, length(circle_list));
            ang_OutDtr      = [];
            for k = 1 : size(transIn_OutDtr, 1)
                pair_k = transIn_OutDtr(k, :);
                if dtrCond == 1
                    dt_Yes = sum((dt_nodes == pair_k(3)), 2);
                else
                    no_control = arrayfun(@(x) ~any(ismember(dt_nodes(x, :), exclude_OutDtr)), (1 : length(dt_nodes))');
                    dt_Yes = sum((dt_nodes == pair_k(3)), 2) & no_control;
                end
                % ------ RTs ------
                RTs_OutDtr = [RTs_OutDtr; respRT_Col(find(from_nodes == pair_k(1) & to_nodes == pair_k(2) & dt_Yes == 1 & errorIdx_Col ~= 1 & fast_Col ~= 1))];
                % ------ choice accuracy ------
                acc_OutDtr = [acc_OutDtr; length(find(choiceId_ang == 1 & from_nodes == pair_k(1) & to_nodes == pair_k(2) & dt_Yes == 1)), length(find(~isnan(choiceId_ang) & from_nodes == pair_k(1) & to_nodes == pair_k(2) & dt_Yes == 1))];
                % ------ choice accuracy trajectory ------
                for iTp = 1 : length(circle_list)
                    choiceId_iTp = choiceId(:, iTp);
                    acc_OutDtr_traj{iTp} = [acc_OutDtr_traj{iTp}; length(find(choiceId_iTp == 1 & from_nodes == pair_k(1) & to_nodes == pair_k(2) & dt_Yes == 1)), length(find(~isnan(choiceId_iTp) & from_nodes == pair_k(1) & to_nodes == pair_k(2) & dt_Yes == 1))];
                end
                % ------ quantify the distance between cue/target and the non-lure distractor ------
                cue_j = pair_k(1);
                tgt_j = pair_k(2);
                dtr_j = dt_nodes(from_nodes == pair_k(1) & to_nodes == pair_k(2) & dt_Yes == 1, :);
                for dd = 1 : size(dtr_j, 1)
                    dtr_j_dd = dtr_j(dd, :);
                    dtr_j_dd = dtr_j_dd(dtr_j_dd ~= 0 & dtr_j_dd ~= pair_k(3));
                    minD_cue_i = miniDmat(dtr_j_dd', cue_j);
                    minD_tgt_i = miniDmat(dtr_j_dd', tgt_j);

                    cue_dtr_distance_subj{SubIdx, 2, iExp} = [cue_dtr_distance_subj{SubIdx, 2, iExp}; minD_cue_i];
                    tgt_dtr_distance_subj{SubIdx, 2, iExp} = [tgt_dtr_distance_subj{SubIdx, 2, iExp}; minD_tgt_i];
                end
                % ------ quantify the relative radian between ang(mouse-target) and ang(mouse-lure/control) ------
                ang_tgt      = angleMag_oneH(find(from_nodes == pair_k(1) & to_nodes == pair_k(2) & dt_Yes == 1), 1);
                ang_lure_tmp = angleMag_oneH(find(from_nodes == pair_k(1) & to_nodes == pair_k(2) & dt_Yes == 1), 2 : end);
                dt_nodes_tmp = dt_nodes(find(from_nodes == pair_k(1) & to_nodes == pair_k(2) & dt_Yes == 1), :);
                ang_lure_tmp = ang_lure_tmp(~isnan(ang_tgt), :);
                dt_nodes_tmp = dt_nodes_tmp(~isnan(ang_tgt), :);
                ang_tgt      = ang_tgt(~isnan(ang_tgt));
                ang_lure     = nan(size(ang_lure_tmp, 1), 1);
                for ll = 1 : size(ang_lure_tmp, 1)
                    ang_lure_tmp_ll = ang_lure_tmp(ll, :);
                    ang_lure(ll) = ang_lure_tmp_ll(find(dt_nodes_tmp(ll, :) == pair_k(3)));
                end
                ang_OutDtr = [ang_OutDtr; ang_tgt, ang_lure];
            end
            if ~isempty(ang_OutDtr)
                ang_OutDtr_diff = mod(ang_OutDtr(:, 2) - ang_OutDtr(:, 1) + pi, 2*pi) - pi;
                ang_transIn_subj(i, 2, SubIdx, iExp) = atan2(mean(sin(ang_OutDtr_diff)), mean(cos(ang_OutDtr_diff)));
            end

            RT_transIn_subj(i, 2, SubIdx, iExp)  = nanmean(RTs_OutDtr);
            acc_transIn_subj(i, 2, SubIdx, iExp) = sum(acc_OutDtr(:, 1)) / sum(acc_OutDtr(:, 2));

            for iTp = 1 : length(circle_list)
                % acc_transIn_traj_subj = zeros(bndNode_Num, length(circle_list), 2, subLen, length(expList));
                % ------ with lure stimulus ------
                acc_transIn_traj_subj(i, iTp, 1, SubIdx, iExp) = sum(acc_InDtr_traj{iTp}(:, 1)) / sum(acc_InDtr_traj{iTp}(:, 2));

                % ------ without lure stimulus
                acc_transIn_traj_subj(i, iTp, 2, SubIdx, iExp) = sum(acc_OutDtr_traj{iTp}(:, 1)) / sum(acc_OutDtr_traj{iTp}(:, 2));
            end

            %%% trial length of each combination
            transIn_Len(i, 1, SubIdx, iExp) = length(RTs_InDtr);
            transIn_Len(i, 2, SubIdx, iExp) = length(RTs_OutDtr);

            %% ---------transition: from boundary node to boundary node---------
            transOut_InDtr  = transOut_clsInDtr{i};
            transOut_OutDtr = transOut_clsOutDtr{i};
            exclude_InDtr   = unique(transOut_OutDtr(:, 3)); % should be 3 nodes
            exclude_OutDtr  = unique(transOut_InDtr(:, 3));  % should be 1 node
            % **************** 3. boundary-to-boundary transition: with lure ****************
            RTs_InDtr = [];
            acc_InDtr = [];
            acc_InDtr_traj = cell(1, length(circle_list));
            ang_InDtr      = [];
            for j = 1 : size(transOut_InDtr, 1)
                pair_j = transOut_InDtr(j, :);
                if dtrCond == 1
                    dt_Yes = sum((dt_nodes == pair_j(3)), 2); 
                else
                    no_control = arrayfun(@(x) ~any(ismember(dt_nodes(x, :), exclude_InDtr)), (1 : length(dt_nodes))');
                    dt_Yes = sum((dt_nodes == pair_j(3)), 2) & no_control;
                end
                % ------ RTs ------
                RTs_InDtr = [RTs_InDtr; respRT_Col(find(from_nodes == pair_j(1) & to_nodes == pair_j(2) & dt_Yes == 1 & errorIdx_Col ~= 1 & fast_Col ~= 1))];
                % ------ choice accuracy ------
                acc_InDtr = [acc_InDtr; length(find(choiceId_ang == 1 & from_nodes == pair_j(1) & to_nodes == pair_j(2) & dt_Yes == 1)), length(find(~isnan(choiceId_ang) & from_nodes == pair_j(1) & to_nodes == pair_j(2) & dt_Yes == 1))];
                % ------ choice accuracy trajectory ------
                for iTp = 1 : length(circle_list)
                    choiceId_iTp = choiceId(:, iTp);
                    acc_InDtr_traj{iTp} = [acc_InDtr_traj{iTp}; length(find(choiceId_iTp == 1 & from_nodes == pair_j(1) & to_nodes == pair_j(2) & dt_Yes == 1)), length(find(~isnan(choiceId_iTp) & from_nodes == pair_j(1) & to_nodes == pair_j(2) & dt_Yes == 1))];
                end
                % ------ quantify the distance between cue/target and the non-lure distractor ------
                cue_j = pair_j(1);
                tgt_j = pair_j(2);
                dtr_j = dt_nodes(from_nodes == pair_j(1) & to_nodes == pair_j(2) & dt_Yes == 1, :);
                for dd = 1 : size(dtr_j, 1)
                    dtr_j_dd = dtr_j(dd, :);
                    dtr_j_dd = dtr_j_dd(dtr_j_dd ~= 0 & dtr_j_dd ~= pair_j(3));
                    minD_cue_i = miniDmat(dtr_j_dd', cue_j);
                    minD_tgt_i = miniDmat(dtr_j_dd', tgt_j);

                    cue_dtr_distance_subj{SubIdx, 3, iExp} = [cue_dtr_distance_subj{SubIdx, 3, iExp}; minD_cue_i];
                    tgt_dtr_distance_subj{SubIdx, 3, iExp} = [tgt_dtr_distance_subj{SubIdx, 3, iExp}; minD_tgt_i];
                end
                % ------ quantify the relative radian between ang(mouse-target) and ang(mouse-lure/control) ------
                ang_tgt      = angleMag_oneH(find(from_nodes == pair_j(1) & to_nodes == pair_j(2) & dt_Yes == 1), 1);
                ang_lure_tmp = angleMag_oneH(find(from_nodes == pair_j(1) & to_nodes == pair_j(2) & dt_Yes == 1), 2 : end);
                dt_nodes_tmp = dt_nodes(find(from_nodes == pair_j(1) & to_nodes == pair_j(2) & dt_Yes == 1), :);
                ang_lure_tmp = ang_lure_tmp(~isnan(ang_tgt), :);
                dt_nodes_tmp = dt_nodes_tmp(~isnan(ang_tgt), :);
                ang_tgt      = ang_tgt(~isnan(ang_tgt));
                ang_lure     = nan(size(ang_lure_tmp, 1), 1);
                for ll = 1 : size(ang_lure_tmp, 1)
                    ang_lure_tmp_ll = ang_lure_tmp(ll, :);
                    ang_lure(ll) = ang_lure_tmp_ll(find(dt_nodes_tmp(ll, :) == pair_j(3)));
                end
                ang_InDtr = [ang_InDtr; ang_tgt, ang_lure];
            end
            if ~isempty(ang_InDtr)
                ang_InDtr_diff = mod(ang_InDtr(:, 2) - ang_InDtr(:, 1) + pi, 2*pi) - pi;
                ang_transOut_subj(i, 1, SubIdx, iExp) = atan2(mean(sin(ang_InDtr_diff)), mean(cos(ang_InDtr_diff)));
            end

            RT_transOut_subj(i, 1, SubIdx, iExp)  = nanmean(RTs_InDtr);
            acc_transOut_subj(i, 1, SubIdx, iExp) = sum(acc_InDtr(:, 1)) / sum(acc_InDtr(:, 2));

            % **************** 4. boundary-to-boundary transition: without lure ****************
            RTs_OutDtr = [];
            acc_OutDtr = [];
            acc_OutDtr_traj = cell(1, length(circle_list));
            ang_OutDtr      = [];
            ang_OutDtr_traj = cell(1, length(circle_list));
            for k = 1 : size(transOut_OutDtr, 1)
                pair_k = transOut_OutDtr(k, :);
                if dtrCond == 1
                    dt_Yes = sum((dt_nodes == pair_k(3)), 2);
                else
                    no_control = arrayfun(@(x) ~any(ismember(dt_nodes(x, :), exclude_OutDtr)), (1 : length(dt_nodes))');
                    dt_Yes = sum((dt_nodes == pair_k(3)), 2) & no_control;
                end
                % ------ RTs ------
                RTs_OutDtr = [RTs_OutDtr; respRT_Col(find(from_nodes == pair_k(1) & to_nodes == pair_k(2) & dt_Yes == 1 & errorIdx_Col ~= 1 & fast_Col ~= 1))];
                % ------ choice accuracy ------
                acc_OutDtr = [acc_OutDtr; length(find(choiceId_ang == 1 & from_nodes == pair_k(1) & to_nodes == pair_k(2) & dt_Yes == 1)), length(find(~isnan(choiceId_ang) & from_nodes == pair_k(1) & to_nodes == pair_k(2) & dt_Yes == 1))];
                % ------ choice accuracy trajectory ------
                for iTp = 1 : length(circle_list)
                    choiceId_iTp = choiceId(:, iTp);
                    acc_OutDtr_traj{iTp} = [acc_OutDtr_traj{iTp}; length(find(choiceId_iTp == 1 & from_nodes == pair_k(1) & to_nodes == pair_k(2) & dt_Yes == 1)), length(find(~isnan(choiceId_iTp) & from_nodes == pair_k(1) & to_nodes == pair_k(2) & dt_Yes == 1))];
                end
                % ------ quantify the distance between cue/target and the non-lure distractor ------
                cue_j = pair_k(1);
                tgt_j = pair_k(2);
                dtr_j = dt_nodes(from_nodes == pair_k(1) & to_nodes == pair_k(2) & dt_Yes == 1, :);
                for dd = 1 : size(dtr_j, 1)
                    dtr_j_dd = dtr_j(dd, :);
                    dtr_j_dd = dtr_j_dd(dtr_j_dd ~= 0 & dtr_j_dd ~= pair_k(3));
                    minD_cue_i = miniDmat(dtr_j_dd', cue_j);
                    minD_tgt_i = miniDmat(dtr_j_dd', tgt_j);

                    cue_dtr_distance_subj{SubIdx, 4, iExp} = [cue_dtr_distance_subj{SubIdx, 4, iExp}; minD_cue_i];
                    tgt_dtr_distance_subj{SubIdx, 4, iExp} = [tgt_dtr_distance_subj{SubIdx, 4, iExp}; minD_tgt_i];
                end
                % ------ quantify the relative radian between ang(mouse-target) and ang(mouse-lure/control) ------
                ang_tgt      = angleMag_oneH(find(from_nodes == pair_k(1) & to_nodes == pair_k(2) & dt_Yes == 1), 1);
                ang_lure_tmp = angleMag_oneH(find(from_nodes == pair_k(1) & to_nodes == pair_k(2) & dt_Yes == 1), 2 : end);
                dt_nodes_tmp = dt_nodes(find(from_nodes == pair_k(1) & to_nodes == pair_k(2) & dt_Yes == 1), :);
                ang_lure_tmp = ang_lure_tmp(~isnan(ang_tgt), :);
                dt_nodes_tmp = dt_nodes_tmp(~isnan(ang_tgt), :);
                ang_tgt      = ang_tgt(~isnan(ang_tgt));
                ang_lure     = nan(size(ang_lure_tmp, 1), 1);
                for ll = 1 : size(ang_lure_tmp, 1)
                    ang_lure_tmp_ll = ang_lure_tmp(ll, :);
                    ang_lure(ll) = ang_lure_tmp_ll(find(dt_nodes_tmp(ll, :) == pair_k(3)));
                end
                ang_OutDtr = [ang_OutDtr; ang_tgt, ang_lure];
            end
            if ~isempty(ang_OutDtr)
                ang_OutDtr_diff = mod(ang_OutDtr(:, 2) - ang_OutDtr(:, 1) + pi, 2*pi) - pi;
                ang_transOut_subj(i, 2, SubIdx, iExp) = atan2(mean(sin(ang_OutDtr_diff)), mean(cos(ang_OutDtr_diff)));
            end

            RT_transOut_subj(i, 2, SubIdx, iExp)  = nanmean(RTs_OutDtr);
            acc_transOut_subj(i, 2, SubIdx, iExp) = sum(acc_OutDtr(:, 1)) / sum(acc_OutDtr(:, 2));

            for iTp = 1 : length(circle_list)
                % acc_transOut_traj_subj = zeros(bndNode_Num, length(circle_list), 2, subLen, length(expList));
                % ------ with lure stimulus ------
                acc_transOut_traj_subj(i, iTp, 1, SubIdx, iExp) = sum(acc_InDtr_traj{iTp}(:, 1)) / sum(acc_InDtr_traj{iTp}(:, 2));

                % ------ without lure stimulus
                acc_transOut_traj_subj(i, iTp, 2, SubIdx, iExp) = sum(acc_OutDtr_traj{iTp}(:, 1)) / sum(acc_OutDtr_traj{iTp}(:, 2));
            end

            %%% trial length of each combination
            transOut_Len(i, 1, SubIdx, iExp) = length(RTs_InDtr);
            transOut_Len(i, 2, SubIdx, iExp) = length(RTs_OutDtr);

        end

        %% ---------- A new analysis for the boundary-to-within or boundary-to-boundary analysis ----------
        % The basic idea is: for the cue node as boundary in a trial, we
        % can check the condition where the lure stimulus and a 2-step
        % control stimulus exist simultaneously, for these trials, if
        % partcipants responded incorrectly whether they hava more
        % probability to select the lure over the 2-step control


    end
end

%% color settings
colorSet = [249, 183, 176; ... % red 
            84, 185, 211; ...  % blue
            138, 170, 51; ...  % green
            248, 218, 172; ...
            184, 204, 225; ...
            210, 234, 200; ...
            198, 127, 192; ...
            219, 204, 226] ./ [255, 255, 255];
redGrad  = [189, 0, 38; ...
            240, 59, 32; ...
            253, 141, 60; ...
            154, 178, 76; ...
            254, 217, 118; ...
            255, 255, 178] ./ 255;
blueGrad = [8, 81, 156; ...
            49, 130, 189; ...
            107, 174, 214; ...
            158, 202, 225; ...
            198, 219, 239; ...
            239, 243, 255] ./ 255;
greeGrad = [0, 104, 55; ...
            49, 163, 84; ...
            120, 198, 121; ...
            173, 221, 142; ...
            217, 240, 163; ...
            255, 255, 204] ./ 255;


%% average
% % ---------- RTs ----------
% RT_transIn_subj  = zeros(bndNode_Num, 2, subLen, length(expList)); % 2: with vs. without lure stimulus
% RT_transOut_subj = zeros(bndNode_Num, 2, subLen, length(expList));
% % ---------- Choice accuracy: one-hot measure ----------
% acc_transIn_subj  = zeros(bndNode_Num, 2, subLen, length(expList));
% acc_transOut_subj = zeros(bndNode_Num, 2, subLen, length(expList));
% % ---------- Choice accuracy: trajectory ----------
% acc_transIn_traj_subj  = zeros(bndNode_Num, length(circle_list), 2, subLen, length(expList));
% acc_transOut_traj_subj = zeros(bndNode_Num, length(circle_list), 2, subLen, length(expList));
% ---------- Relative angles between the mouse and the target/lure-distractor/control ----------
% ang_transIn_subj       = zeros(bndNode_Num, 2, subLen, length(expList));
% ang_transOut_subj      = zeros(bndNode_Num, 2, subLen, length(expList));

% ---------- average across the boundar nodes ----------
% ------ RTs ------
RT_transIn_nodeAvg     = squeeze(nanmean(RT_transIn_subj, 1)); % 2 * subLen * length(expList)
RT_transOut_nodeAvg    = squeeze(nanmean(RT_transOut_subj, 1));
RT_transInDif_nodeAvg  = squeeze(nanmean(RT_transIn_subj(:, 2, :) - RT_transIn_subj(:, 1, :), 1)); % 1 * subLen
RT_transOutDif_nodeAvg = squeeze(nanmean(RT_transOut_subj(:, 2, :) - RT_transOut_subj(:, 1, :), 1)); 
% ------ choice accuracy ------
acc_transIn_nodeAvg    = squeeze(nanmean(acc_transIn_subj, 1));
acc_transOut_nodeAvg   = squeeze(nanmean(acc_transOut_subj, 1));
% ------ choice accuracy trajectory ------
acc_transIn_traj_nodeAvg  = squeeze(nanmean(acc_transIn_traj_subj, 1)); % length(circle_list) * 2 * subLen * length(expList)) 
acc_transOut_traj_nodeAvg = squeeze(nanmean(acc_transOut_traj_subj, 1));
% ------ relative angles ------
ang_transIn_nodeAvg  = squeeze(...
    atan2(mean(sin(ang_transIn_subj), 1), ...
          mean(cos(ang_transIn_subj), 1)));
ang_transOut_nodeAvg = squeeze(...
    atan2(mean(sin(ang_transOut_subj), 1), ...
          mean(cos(ang_transOut_subj), 1)));

%% which data to plot
dataFlg = 2;
if dataFlg == 1    % RTs
    trans_boundary_within   = RT_transIn_nodeAvg;
    trans_boundary_boundary = RT_transOut_nodeAvg;
elseif dataFlg == 2 % choice accuracy
    trans_boundary_within   = acc_transIn_nodeAvg;
    trans_boundary_boundary = acc_transOut_nodeAvg;
end

%% BehavioralPaper, Figure xx: influence of lure stimulus on choice accuracy (one-hot metric)
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
tail_id  = 'left'; %'both';
if dtrCond == 1
    chance_i = 1/2;
elseif dtrCond == 2
    chance_i = 1/3;
elseif dtrCond == 3
    chance_i = 1/4;
elseif dtrCond ==4
    chance_i = mean([1/2, 1/3, 1/4]);
end
for iExp = 1 : 3
    disp(['---------- ', expList{iExp}, ' ----------']);
    figure('Position', [100 100 260 120]), clf;

    barPos = [1, 1.5; 1.7, 2.2];
    for iTb = 1 : 2     % 'transition from boundary node to within node vs. from boundary to boundary'
        if iTb == 1
            transWord     = 'boundary-to-within';
            transData_iTb = trans_boundary_within(:, :, iExp);
        elseif iTb == 2
            transWord     = 'boundary-to-boundary';
            transData_iTb = trans_boundary_boundary(:, :, iExp);
        end
        transData_iTb = transData_iTb'; % subLen * 2 (with vs. without lure stimulus)
        [accAvg_rOh, accSem_rOh] = Mean_and_Se(transData_iTb, 1);

        barPos_i = barPos(iTb, :);
        %%% line plot
        plot(barPos_i, transData_iTb, 'Color', [0.6, 0.6, 0.6], 'LineStyle', '-', 'LineWidth', indvLineW); hold on;
        plot(barPos_i, accAvg_rOh, 'Color', [0, 0, 0], 'LineStyle', '-', 'LineWidth', errLineWid); hold on;
        for ilr = 1 : 2 % 'lure distractor exists vs. none'
            if ilr == 1
                colorTmp = [0, 0, 0]; % with lure distractor
            elseif ilr == 2
                colorTmp = [1, 1, 1]; % without lure distractor
            end
            errorbar(barPos_i(ilr), accAvg_rOh(ilr), accSem_rOh(ilr), 'Color', 'k', 'LineStyle', 'none', 'LineWidth', errLineWid); hold on;
            plot(barPos_i(ilr), accAvg_rOh(ilr), 'Marker', 'o', 'MarkerSize', markSize, 'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', colorTmp, 'LineStyle', '-'); hold on;
        end
        % ------ Statistical tests ------
        disp(['======== ', transWord, ': with vs. without lure stimulus ========']);
        [h, p, ci, stats] = ttest(transData_iTb(:, 1), transData_iTb(:, 2), 'tail', tail_id)
        disp(['t=', num2str(stats.tstat, '%4.3f'), ', p=', num2str(p, '%4.3f')])
    end
    xlim([0.6, 2.6]);
    if dataFlg == 1    % RTs
        ylim([0.5, 2]);
        if figKey == 0
            % ------For presentation------
            set(gca, 'LineWidth', 2);
            set(gca, 'FontSize', 15, 'FontWeight', 'bold', 'FontName', 'Arial');
            set(gca, 'XTick', '', 'XTickLabel', '');
            set(gca, 'YTick', [0.5, 1, 2], 'YTickLabel', [0.5, 1, 2]);
        elseif figKey == 1
            % ------For Adobe Illustrator------
            set(gca, 'LineWidth', 0.8);
            set(gca, 'FontSize', 10, 'FontWeight', 'bold', 'FontName', 'Arial');
            set(gca, 'XTick', [1, 1.5, 1.7, 2.2], 'XTickLabel', '');
            set(gca, 'YTick', [0.5, 1, 2], 'YTickLabel', '');
        end

    elseif dataFlg == 2 % choice accuracy
        ylim([0, 1]);
        %ylim([0.05, 0.8]);
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
            set(gca, 'YTick', [0.05, 0.2 : 0.2 : 0.8], 'YTickLabel', '');
        end
    end
   
    box off;
    ax = gca;
    save_name = ['Exp', num2str(iExp), '-', expList{iExp}, '-oneHot-lure.png'];
%     exportgraphics(ax, save_name, 'Resolution', 600);
end

%% BehavioralPaper, Figure xx: influence of lure stimulus on choice accuracy trajectory
% acc_transIn_traj_nodeAvg  = squeeze(nanmean(acc_transIn_traj_subj, 1)); % length(circle_list) * 2 * subLen * length(expList)) 
% acc_transOut_traj_nodeAvg = squeeze(nanmean(acc_transOut_traj_subj, 1));
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
if dtrCond == 1
    chance_i = 1/2;
elseif dtrCond == 2
    chance_i = 1/3;
elseif dtrCond == 3
    chance_i = 1/4;
elseif dtrCond ==4
    chance_i = mean([1/2, 1/3, 1/4]);
end

for iExp = 1 : 3
    disp(['---------- ', expList{iExp}, ' ----------']);
    figure('Position', [100 100 260 120]), clf;

    for iTb = 1 : 2     % 'transition from boundary node to within node vs. from boundary to boundary'
        if iTb == 1
            transWord          = 'boundary-to-within';
            transData_traj_iTb = acc_transIn_traj_nodeAvg(:, :, :, iExp); % length(circle_list) * 2 * subLen
            color_iTb          = redGrad(2, :);
        elseif iTb == 2
            transWord          = 'boundary-to-boundary';
            transData_traj_iTb = acc_transOut_traj_nodeAvg(:, :, :, iExp);
            color_iTb          = blueGrad(2, :);
        end
        [accAvg_rOh, accSem_rOh] = Mean_and_Se(transData_traj_iTb, 3); % length(circle_list) * 2
        for ilr = 1 : 2 % 'lure distractor exists vs. none'
            if ilr == 1
                LineSty = '-';
            elseif ilr == 2
                LineSty = ':';
            end
            shadedErrorBar(circle_list, accAvg_rOh(:, ilr), accSem_rOh(:, ilr), {'Color', color_iTb, 'MarkerFaceColor', color_iTb, 'LineStyle', LineSty, 'LineWidth', 2}, 0.5); hold on;
        end
    end
    xlim([0, 0.8]);
    ylim([0, 0.6]);
    plot(xlim, [chance_i, chance_i], 'k--', 'LineWidth', 0.8); hold on;
    if iExp == 1 || iExp == 3
        plot([0.8, 0.8], ylim, 'k--', 'LineWidth', 1); hold on;
    end
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
        set(gca, 'XTick', [0, 0.4, 0.8], 'XTickLabel', '');
        set(gca, 'YTick', [0, 0.3, 0.6], 'YTickLabel', '');
    end
    box off;
    ax = gca;
    save_name = ['Exp', num2str(iExp), '-', expList{iExp}, '-learning-lure.png'];
    %exportgraphics(ax, save_name, 'Resolution', 600);
end

%% quantify the distance between the cue and non-lure/control stimulus or distance between the target and non-lure/control stimulus
% cue_dtr_distance_subj = cell(subLen, 4, length(expList)); % 4: with vs. without lure stimulus for boundary-to-within and boundary-to-boundary transitions
% tgt_dtr_distance_subj = cell(subLen, 4, length(expList));
cue_dtr_disAvg = nan(subLen, 4, length(expList));
tgt_dtr_disAvg = nan(subLen, 4, length(expList));
for iExp = 1 : 3
    for SubIdx = 1 : subLen
        for ii = 1 : 4
            % ------ distance between cue and non-lure/control stimulus ------
            cue_dtr_disAvg(SubIdx, ii, iExp) = nanmean(cue_dtr_distance_subj{SubIdx, ii, iExp});

            % ------ distance between target and non-lure/control stimulus ------
            tgt_dtr_disAvg(SubIdx, ii, iExp) = nanmean(tgt_dtr_distance_subj{SubIdx, ii, iExp});
        end
    end
end
% ------ which data to plot ------
disFlg = 1;
if disFlg == 1     % distance between cue and non-lure stimulus
    distance_cal = cue_dtr_disAvg;
elseif disFlg == 2 % distance between target and non-lure stimulus
    distance_cal = tgt_dtr_disAvg;
end


%% plotting
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
for iExp = 1 : 3
    disp(['---------- ', expList{iExp}, ' ----------']);
    figure('Position', [100 100 260 120]), clf;

    barPos = [1, 1.5; 1.7, 2.2];
    for iTb = 1 : 2     % 'transition from boundary node to within node vs. from boundary to boundary'
        if iTb == 1
            transWord = 'boundary-to-within';
            dis_iTb   = distance_cal(:, 1 : 2, iExp);
        elseif iTb == 2
            transWord = 'boundary-to-boundary';
            dis_iTb   = distance_cal(:, 3 : 4, iExp);
        end
        [disAvg_rOh, disSem_rOh] = Mean_and_Se(dis_iTb, 1);

        barPos_i = barPos(iTb, :);
        %%% line plot
        plot(barPos_i, dis_iTb, 'Color', [0.6, 0.6, 0.6], 'LineStyle', '-', 'LineWidth', indvLineW); hold on;
        plot(barPos_i, disAvg_rOh, 'Color', [0, 0, 0], 'LineStyle', '-', 'LineWidth', errLineWid); hold on;
        for ilr = 1 : 2 % 'lure distractor exists vs. none'
            if ilr == 1
                colorTmp = [0, 0, 0]; % with lure distractor
            elseif ilr == 2
                colorTmp = [1, 1, 1]; % without lure distractor
            end
            errorbar(barPos_i(ilr), disAvg_rOh(ilr), disSem_rOh(ilr), 'Color', 'k', 'LineStyle', 'none', 'LineWidth', errLineWid); hold on;
            plot(barPos_i(ilr), disAvg_rOh(ilr), 'Marker', 'o', 'MarkerSize', markSize, 'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', colorTmp, 'LineStyle', '-'); hold on;
        end
        % ------ Statistical tests ------
        disp(['======== ', transWord, ': with vs. without lure stimulus ========']);
        [h, p, ci, stats] = ttest(dis_iTb(:, 1), dis_iTb(:, 2))
        disp(['t=', num2str(stats.tstat, '%4.3f'), ', p=', num2str(p, '%4.3f')])
    end
    xlim([0.6, 2.6]);
    if dtrCond ==4
        %ylim([2.5, 3.5]);
    end
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
        set(gca, 'YTick', [2.5, 3, 3.5], 'YTickLabel', '');
    end
    box off;
    ax = gca;
    save_name = ['Exp', num2str(iExp), '-', expList{iExp}, '-distance-lure.png'];
%     exportgraphics(ax, save_name, 'Resolution', 600);
end

%% The relative angle difference between mouse-target and mouse-lure/control
% ---------- Relative angles between the mouse and the target/lure-distractor/control ----------
% ang_transIn_subj     = zeros(bndNode_Num, 2, subLen, length(expList));
% ang_transOut_subj    = zeros(bndNode_Num, 2, subLen, length(expList));
% ------ relative angles ------
% ang_transIn_nodeAvg  = squeeze(...
%     atan2(mean(sin(ang_transIn_subj), 1), ...
%           mean(cos(ang_transIn_subj), 1)));
% ang_transOut_nodeAvg = squeeze(...
%     atan2(mean(sin(ang_transOut_subj), 1), ...
%           mean(cos(ang_transOut_subj), 1)));

% figKey = 1;
% if figKey == 0
%     barLineWid = 2;
%     errLineWid = 3;
%     refLineWid = 1;
%     indvLineW  = 1;
%     markSize   = 6;
% elseif figKey == 1
%     barLineWid = 1;
%     errLineWid = 2;
%     refLineWid = 0.5;
%     indvLineW  = 0.4;
%     markSize   = 4.5;
% end
% for iExp = 1 : 3
%     disp(['---------- ', expList{iExp}, ' ----------']);
%     figure('Position', [100 100 260 120]), clf;
% 
%     barPos = [1, 1.5; 1.7, 2.2];
%     for iTb = 1 : 2     % 'transition from boundary node to within node vs. from boundary to boundary'
%         if iTb == 1
%             transWord     = 'boundary-to-within';
%             transData_iTb = ang_transIn_nodeAvg(:, :, iExp);
%         elseif iTb == 2
%             transWord     = 'boundary-to-boundary';
%             transData_iTb = ang_transOut_nodeAvg(:, :, iExp);
%         end
%         transData_iTb = transData_iTb'; % subLen * 2 (with vs. without lure stimulus)
%         % Number of valid samples per condition
%         N = sum(~isnan(transData_iTb), 1);   % 1 × 2
% 
%         % Circular mean (NaN-safe)
%         sinMean = mean(sin(transData_iTb), 1, 'omitnan');
%         cosMean = mean(cos(transData_iTb), 1, 'omitnan');
% 
%         circMean_sub = atan2(sinMean, cosMean);
% 
%         % Resultant vector length
%         R = sqrt(sinMean.^2 + cosMean.^2);
% 
%         % Circular SEM (uses effective N)
%         circSEM = sqrt(-2 * log(R)) ./ sqrt(N);
% 
%         circMean_sub = squeeze(circMean_sub);   % 1 × 2
%         circSEM      = squeeze(circSEM);        % 1 × 2
% 
%         barPos_i = barPos(iTb, :);
%         %%% line plot
%         plot(barPos_i, transData_iTb, 'Color', [0.6, 0.6, 0.6], 'LineStyle', '-', 'LineWidth', indvLineW); hold on;
%         plot(barPos_i, circMean_sub, 'Color', [0, 0, 0], 'LineStyle', '-', 'LineWidth', errLineWid); hold on;
%         for ilr = 1 : 2 % 'lure distractor exists vs. none'
%             if ilr == 1
%                 colorTmp = [0, 0, 0]; % with lure distractor
%             elseif ilr == 2
%                 colorTmp = [1, 1, 1]; % without lure distractor
%             end
%             errorbar(barPos_i(ilr), circMean_sub(ilr), circSEM(ilr), 'Color', 'k', 'LineStyle', 'none', 'LineWidth', errLineWid); hold on;
%             plot(barPos_i(ilr), circMean_sub(ilr), 'Marker', 'o', 'MarkerSize', markSize, 'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', colorTmp, 'LineStyle', '-'); hold on;
%         end
%         % ------ Statistical tests ------
%         disp(['======== ', transWord, ': with vs. without lure stimulus ========']);
%         [h, p, ci, stats] = ttest(transData_iTb(:, 1), transData_iTb(:, 2))
%         disp(['t=', num2str(stats.tstat, '%4.3f'), ', p=', num2str(p, '%4.3f')])
%     end
%     xlim([0.6, 2.6]);
%     %ylim([0, 1]);
%     if figKey == 0
%         % ------For presentation------
%         set(gca, 'LineWidth', 2);
%         set(gca, 'FontSize', 15, 'FontWeight', 'bold', 'FontName', 'Arial');
%         set(gca, 'XTick', '', 'XTickLabel', '');
%         set(gca, 'YTick', 0 : 0.5 : 1, 'YTickLabel', 0 : 0.5 : 1);
%     elseif figKey == 1
%         % ------For Adobe Illustrator------
%         set(gca, 'LineWidth', 0.8);
%         set(gca, 'FontSize', 10, 'FontWeight', 'bold', 'FontName', 'Arial');
%         set(gca, 'XTick', [1, 1.5, 1.7, 2.2], 'XTickLabel', '');
%         %set(gca, 'YTick', [0, 0.5, 1], 'YTickLabel', '');
%     end
%    
%     box off;
% end



