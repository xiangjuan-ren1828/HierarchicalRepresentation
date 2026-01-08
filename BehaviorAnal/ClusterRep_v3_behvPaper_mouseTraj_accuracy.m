% Original script: 
% ClusterRep_ChoiceModel_accuracy.m under 'clusterRep_modeling' folder
% write by rxj @ 08/26/2021
%
% Calculate the accuracy based on classification of mouse trajectory angles

clear
clc

%%
addpath('tight_subplot/');
addpath(genpath('HierarchicalCluster/'));
addpath('fdr_bh');

% addpath(genpath('/home/ubuntu/Xiangjuan/MyExperiment/HierarchicalCluster/'));

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

subLen = 24;
%%% merge all trials
angAcc_exp          = nan(subLen, length(circle_list), length(expList)); %% accuracy in each time point
lenData_exp         = nan(subLen, length(circle_list), length(expList)); %% data length in each time point
angAcc_dtr_exp      = nan(subLen, length(circle_list), 3, length(expList));
%%% choice accuracy trajectory for within- and between-transitions
angAcc_traj_trans_exp = nan(subLen, 2, length(circle_list), length(expList)); % 2: within and between-transitions
%%% choice accuracy trajectory for within- and between-transitions with
%%% final correct or incorrect response separation
angAcc_traj_trans_cor_exp = nan(subLen, 4, length(circle_list), length(expList)); % Col1-2: final response correct vs. incorrect for the within-transitions; Col2: ~ for between-transitions
dtrNo_traj_trans_cor_exp  = nan(subLen, 4, length(circle_list), length(expList)); % Quantifying the distractor number for each condition
%%% only for the Random trials
angAcc_Rand_exp     = nan(subLen, length(circle_list), length(expList)); %% accuracy in each time point
lenData_Rand_exp    = nan(subLen, length(circle_list), length(expList)); %% data length in each time point
angAcc_dtr_Rand_exp = nan(subLen, length(circle_list), 3, length(expList));
%%% only for the Hamiltonian trials
angAcc_Hami_exp     = nan(subLen, length(circle_list), length(expList)); %% accuracy in each time point
lenData_Hami_exp    = nan(subLen, length(circle_list), length(expList)); %% data length in each time point
angAcc_dtr_Hami_exp = nan(subLen, length(circle_list), 3, length(expList));
%%% accuracy for the within and between transitions for Random and Hamiltonian Walk separately
acc_trans_walk      = nan(subLen, 2, 2, length(expList)); % the 1st and 2nd '2' denote 'within vs. between cluster transition' and 'random and hamiltonian' walk
%%% trial proportions of within-cluster and between-cluster transition under Random and Hamiltonian Walk
trlPro_trans_walk   = nan(subLen, 2, 2, length(expList)); 
%%% accuracy for boundary-to-within and boundary-to-boundary transitions
%%% with versus without lure stimulus
acc_trans_lure      = nan(subLen, 2, 2, length(expList)); % the 1st and 2nd '2' denote 'lure distractor exists vs. none' and 'transition from boundary node to within node vs. from boundary to boundary'
acc_trans_traj_lure = nan(subLen, length(circle_list), 2, 2, length(expList));
% first 2:  with or without lure stimulus
% second 2: boundary-to-within transition or bounary-to-boundary transition
% third 3:  random or hamiltonian walk
acc_trans_walk_lure      = nan(subLen, 2, 2, 2, length(expList));
acc_trans_walk_traj_lure = nan(subLen, length(circle_list), 2, 2, 2, length(expList));
%%% learning curves
% all trials together
angAcc_ln_exp            = nan(subLen, BinL, length(expList));
% different distractor numbers
angAcc_ln_dtr_exp        = nan(subLen, BinL, 3, length(expList));
% within- vs. between-trans
angAcc_trans_ln_exp      = nan(subLen, BinL, 2, length(expList));    % 2: within- vs. between transition
% within- vs. between-trans for Random and Hamiltonian Walk
angAcc_trans_walk_ln_exp = nan(subLen, BinL, 2, 2, length(expList)); % first 2: within- vs. between transition; 2nd 2: random vs. hamiltonian walk

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
    
    %% loop over all trials
    %%% merge all trials
    angAcc          = zeros(subLen, length(circle_list)); %% accuracy in each time point
    lenData         = zeros(subLen, length(circle_list)); %% data length in each time point
    angAcc_dtr      = zeros(subLen, length(circle_list), 3); % 3 means distractor number categories
    %%% only for the Random trials
    angAcc_Rand     = zeros(subLen, length(circle_list)); %% accuracy in each time point
    lenData_Rand    = zeros(subLen, length(circle_list)); %% data length in each time point
    angAcc_dtr_Rand = zeros(subLen, length(circle_list), 3); % 3 means distractor number categories
    %%% only for the Hamiltonian trials
    angAcc_Hami     = zeros(subLen, length(circle_list)); %% accuracy in each time point
    lenData_Hami    = zeros(subLen, length(circle_list)); %% data length in each time point
    angAcc_dtr_Hami = zeros(subLen, length(circle_list), 3); % 3 means distractor number categories
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
        transStyle = (transIn(:, 1) == transIn(:, 2));
        nodesLabel = zeros(nTrials, 1);
        nodesLabel(from_nodes == 1 | from_nodes == 5 | from_nodes == 6 | from_nodes == 10 | from_nodes == 11 | from_nodes == 15) = 1;
        
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
        iCount   = 1;
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
        %% merge all trials
        for iTp = 1 : length(circle_list)
            choiceId_i      = choiceId(:, iTp);
            dtrNoCnt_iSub_i = dtrNoCnt_iSub;
            dtrNoCnt_iSub_i(isnan(choiceId_i)) = [];
            choiceId_i(isnan(choiceId_i))      = [];
            angAcc(iSub, iTp)  = length(find(choiceId_i == 1)) / length(choiceId_i);
            lenData(iSub, iTp) = length(choiceId_i);
            
            %%% separately calculating the accuracy for different distractor
            %%% numbers
            for iDtr = 1 : 3
                angAcc_dtr(iSub, iTp, iDtr) = length(find(choiceId_i == 1 & dtrNoCnt_iSub_i == iDtr)) / length(find(dtrNoCnt_iSub_i == iDtr));
            end
        end
        
        %% only for the Random trials
        used_Trl = 1 : rndTrial;
        for iTp = 1 : length(circle_list)
            choiceId_i      = choiceId(used_Trl, iTp);
            dtrNoCnt_iSub_i = dtrNoCnt_iSub(used_Trl);
            dtrNoCnt_iSub_i(isnan(choiceId_i)) = [];
            choiceId_i(isnan(choiceId_i))     = [];
            angAcc_Rand(iSub, iTp)  = length(find(choiceId_i == 1)) / length(choiceId_i);
            lenData_Rand(iSub, iTp) = length(choiceId_i);

            for iDtr = 1 : 3
                angAcc_dtr_Rand(iSub, iTp, iDtr) = length(find(choiceId_i == 1 & dtrNoCnt_iSub_i == iDtr)) / length(find(dtrNoCnt_iSub_i == iDtr));
            end
        end
        
        %% only for the Hamiltonian trials
        used_Trl = 701 : nTrials;
        for iTp = 1 : length(circle_list)
            choiceId_i      = choiceId(used_Trl, iTp);
            dtrNoCnt_iSub_i = dtrNoCnt_iSub(used_Trl);
            dtrNoCnt_iSub_i(isnan(choiceId_i)) = [];
            choiceId_i(isnan(choiceId_i))     = [];
            angAcc_Hami(iSub, iTp)  = length(find(choiceId_i == 1)) / length(choiceId_i);
            lenData_Hami(iSub, iTp) = length(choiceId_i);

            for iDtr = 1 : 3
                angAcc_dtr_Hami(iSub, iTp, iDtr) = length(find(choiceId_i == 1 & dtrNoCnt_iSub_i == iDtr)) / length(find(dtrNoCnt_iSub_i == iDtr));
            end
        end

        %% trajectory when lure stimulus exists
        % choiceId = nan(length(stim), length(circle_list));
        for iTp = 1 : length(circle_list)
            choiceId_i = choiceId(:, iTp);

            transStyle_iTp = transStyle;
            rndHam_Col_iTp = rndHam_Col;
            nodesLabel_iTp = nodesLabel;
            lureIn_iTp     = lureIn;
            transStyle_iTp(isnan(choiceId_i)) = [];
            rndHam_Col_iTp(isnan(choiceId_i)) = [];
            nodesLabel_iTp(isnan(choiceId_i)) = [];
            lureIn_iTp(isnan(choiceId_i))     = [];
            choiceId_i(isnan(choiceId_i))     = [];

            % -----------------------------------------------------------------------------------------
            % the 1st and 2nd '2' denote 'lure distractor exists vs. none' and 'transition from boundary node to within node vs. from boundary to boundary'
            % acc_trans_traj_lure = zeros(subLen, length(circle_list), 2, 2, length(expList));
            acc_trans_traj_lure(iSub, iTp, 1, 1, iExp) = length(find(choiceId_i == 1 & nodesLabel_iTp == 1 & lureIn_iTp == 1 & transStyle_iTp == 1)) / length(find(nodesLabel_iTp == 1 & lureIn_iTp == 1 & transStyle_iTp == 1));
            % ------ lure stimulus exists in boundary-to-boundary trans ------
            acc_trans_traj_lure(iSub, iTp, 1, 2, iExp) = length(find(choiceId_i == 1 & nodesLabel_iTp == 1 & lureIn_iTp == 1 & transStyle_iTp ~= 1)) / length(find(nodesLabel_iTp == 1 & lureIn_iTp == 1 & transStyle_iTp ~= 1));
            % ------ no lure stimulus in boundary-to-within trans ------
            acc_trans_traj_lure(iSub, iTp, 2, 1, iExp) = length(find(choiceId_i == 1 & nodesLabel_iTp == 1 & lureIn_iTp == 0 & transStyle_iTp == 1)) / length(find(nodesLabel_iTp == 1 & lureIn_iTp == 0 & transStyle_iTp == 1));
            % ------ no lure stimulus in boundary-to-boundary trans ------
            acc_trans_traj_lure(iSub, iTp, 2, 2, iExp) = length(find(choiceId_i == 1 & nodesLabel_iTp == 1 & lureIn_iTp == 0 & transStyle_iTp ~= 1)) / length(find(nodesLabel_iTp == 1 & lureIn_iTp == 0 & transStyle_iTp ~= 1));

            % -----------------------------------------------------------------------------------------
            % first 2:  with or without lure stimulus
            % second 2: boundary-to-within transition or bounary-to-boundary transition
            % third 3:  random or hamiltonian walk
            % acc_trans_walk_traj_lure = zeros(subLen, length(circle_list), 2, 2, 2, length(expList));
            % ------ lure stimulus exists ------
            % boundary-to-within trans, Random walk
            acc_trans_walk_traj_lure(iSub, iTp, 1, 1, 1, iExp) = length(find(choiceId_i == 1 & nodesLabel_iTp == 1 & lureIn_iTp == 1 & transStyle_iTp == 1 & rndHam_Col_iTp == 1)) / length(find(nodesLabel_iTp == 1 & lureIn_iTp == 1 & transStyle_iTp == 1 & rndHam_Col_iTp == 1));
            % boundary-to-boundary trans, Random walk
            acc_trans_walk_traj_lure(iSub, iTp, 1, 2, 1, iExp) = length(find(choiceId_i == 1 & nodesLabel_iTp == 1 & lureIn_iTp == 1 & transStyle_iTp ~= 1 & rndHam_Col_iTp == 1)) / length(find(nodesLabel_iTp == 1 & lureIn_iTp == 1 & transStyle_iTp ~= 1 & rndHam_Col_iTp == 1));
            % boundary-to-within trans, Hamiltonian walk
            acc_trans_walk_traj_lure(iSub, iTp, 1, 1, 2, iExp) = length(find(choiceId_i == 1 & nodesLabel_iTp == 1 & lureIn_iTp == 1 & transStyle_iTp == 1 & rndHam_Col_iTp == 2)) / length(find(nodesLabel_iTp == 1 & lureIn_iTp == 1 & transStyle_iTp == 1 & rndHam_Col_iTp == 2));
            % boundary-to-boundary trans, Hamiltonian walk
            acc_trans_walk_traj_lure(iSub, iTp, 1, 2, 2, iExp) = length(find(choiceId_i == 1 & nodesLabel_iTp == 1 & lureIn_iTp == 1 & transStyle_iTp ~= 1 & rndHam_Col_iTp == 2)) / length(find(nodesLabel_iTp == 1 & lureIn_iTp == 1 & transStyle_iTp ~= 1 & rndHam_Col_iTp == 2));
            % ------ no lure stimulus ------
            acc_trans_walk_traj_lure(iSub, iTp, 2, 1, 1, iExp) = length(find(choiceId_i == 1 & nodesLabel_iTp == 1 & lureIn_iTp == 0 & transStyle_iTp == 1 & rndHam_Col_iTp == 1)) / length(find(nodesLabel_iTp == 1 & lureIn_iTp == 0 & transStyle_iTp == 1 & rndHam_Col_iTp == 1));
            acc_trans_walk_traj_lure(iSub, iTp, 2, 2, 1, iExp) = length(find(choiceId_i == 1 & nodesLabel_iTp == 1 & lureIn_iTp == 0 & transStyle_iTp ~= 1 & rndHam_Col_iTp == 1)) / length(find(nodesLabel_iTp == 1 & lureIn_iTp == 0 & transStyle_iTp ~= 1 & rndHam_Col_iTp == 1));
            acc_trans_walk_traj_lure(iSub, iTp, 2, 1, 2, iExp) = length(find(choiceId_i == 1 & nodesLabel_iTp == 1 & lureIn_iTp == 0 & transStyle_iTp == 1 & rndHam_Col_iTp == 2)) / length(find(nodesLabel_iTp == 1 & lureIn_iTp == 0 & transStyle_iTp == 1 & rndHam_Col_iTp == 2));
            acc_trans_walk_traj_lure(iSub, iTp, 2, 2, 2, iExp) = length(find(choiceId_i == 1 & nodesLabel_iTp == 1 & lureIn_iTp == 0 & transStyle_iTp ~= 1 & rndHam_Col_iTp == 2)) / length(find(nodesLabel_iTp == 1 & lureIn_iTp == 0 & transStyle_iTp ~= 1 & rndHam_Col_iTp == 2));
        end

        %% separating the trajectory based on within- vs. between-transitions & final response correct vs. incorrect
        %%% choice accuracy trajectory for within- and between-transitions
        % angAcc_traj_trans_exp = nan(subLen, 2, length(circle_list), length(expList)); % 2: within and between-transitions
        %%% choice accuracy trajectory for within- and between-transitions with
        %%% final correct or incorrect response separation
        % angAcc_traj_trans_cor_exp = nan(subLen, 4, length(circle_list), length(expList)); % Col1-2: final response correct vs. incorrect for the within-transitions; Col2: ~ for between-transitions
        % dtrNo_traj_trans_cor_exp  = nan(subLen, 4, length(circle_list), length(expList)); % Quantifying the distractor number for each condition

        for iTp = 1 : length(circle_list)
            choiceId_i      = choiceId(:, iTp);
            choiceFinal_iTp = choiceId_ang;
            transStyle_iTp  = transStyle;
            dtrNo_iTp       = dtrNoCnt_iSub;

            % ---------- within vs. between transitions ----------
            del_trial       = isnan(choiceId_i);
            transStyle_iiTp = transStyle_iTp(~del_trial);
            choiceid_iiTp   = choiceId_i(~del_trial);
            angAcc_traj_trans_exp(iSub, 1, iTp, iExp) = length(find(choiceid_iiTp == 1 & transStyle_iiTp == 1)) / length(find(transStyle_iiTp == 1));
            angAcc_traj_trans_exp(iSub, 2, iTp, iExp) = length(find(choiceid_iiTp == 1 & transStyle_iiTp ~= 1)) / length(find(transStyle_iiTp ~= 1));

            % ---------- within vs. between transitions for final response correct vs. incorrect ----------
            del_trial        = isnan(choiceId_i) | isnan(choiceFinal_iTp);
            transStyle_iiTp  = transStyle_iTp(~del_trial);
            choiceid_iiTp    = choiceId_i(~del_trial);
            choiceFinal_iiTp = choiceFinal_iTp(~del_trial);
            angAcc_traj_trans_cor_exp(iSub, 1, iTp, iExp) = length(find(choiceid_iiTp == 1 & transStyle_iiTp == 1 & choiceFinal_iiTp == 1)) / length(find(transStyle_iiTp == 1 & choiceFinal_iiTp == 1));
            angAcc_traj_trans_cor_exp(iSub, 2, iTp, iExp) = length(find(choiceid_iiTp == 1 & transStyle_iiTp == 1 & choiceFinal_iiTp == 0)) / length(find(transStyle_iiTp == 1 & choiceFinal_iiTp == 0));
            angAcc_traj_trans_cor_exp(iSub, 3, iTp, iExp) = length(find(choiceid_iiTp == 1 & transStyle_iiTp ~= 1 & choiceFinal_iiTp == 1)) / length(find(transStyle_iiTp ~= 1 & choiceFinal_iiTp == 1));
            angAcc_traj_trans_cor_exp(iSub, 4, iTp, iExp) = length(find(choiceid_iiTp == 1 & transStyle_iiTp ~= 1 & choiceFinal_iiTp == 0)) / length(find(transStyle_iiTp ~= 1 & choiceFinal_iiTp == 0));

            % ---------- For each within vs. between & final response
            % correct vs. incorrect condition, quantifying the distractor
            % numbers ----------
            dtrNo_iiTp = dtrNo_iTp(~del_trial);
            if ~isempty(find(choiceid_iiTp == 1 & transStyle_iiTp == 1 & choiceFinal_iiTp == 1))
                dtrNo_traj_trans_cor_exp(iSub, 1, iTp, iExp) = nanmean(dtrNo_iiTp(find(choiceid_iiTp == 1 & transStyle_iiTp == 1 & choiceFinal_iiTp == 1)) + 1);
            end
            if ~isempty(find(choiceid_iiTp == 1 & transStyle_iiTp == 1 & choiceFinal_iiTp == 0))
                dtrNo_traj_trans_cor_exp(iSub, 2, iTp, iExp) = nanmean(dtrNo_iiTp(find(choiceid_iiTp == 1 & transStyle_iiTp == 1 & choiceFinal_iiTp == 0)) + 1);
            end
            if ~isempty(find(choiceid_iiTp == 1 & transStyle_iiTp ~= 1 & choiceFinal_iiTp == 1))
                dtrNo_traj_trans_cor_exp(iSub, 3, iTp, iExp) = nanmean(dtrNo_iiTp(find(choiceid_iiTp == 1 & transStyle_iiTp ~= 1 & choiceFinal_iiTp == 1)) + 1);
            end
            if ~isempty(find(choiceid_iiTp == 1 & transStyle_iiTp ~= 1 & choiceFinal_iiTp == 0))
                dtrNo_traj_trans_cor_exp(iSub, 4, iTp, iExp) = nanmean(dtrNo_iiTp(find(choiceid_iiTp == 1 & transStyle_iiTp ~= 1 & choiceFinal_iiTp == 0)) + 1);
            end

        end

        %% learning effect: accuracy in sliding bins at all sampling time points
        % ------ all trials together ------
        % angAcc_ln_exp            = nan(subLen, BinL, length(expList));
        for iB = 1 : BinL
            trlIdx = (iB - 1) * nBin + 1 : iB * nBin;
            choiceId_trl = choiceId_ang(trlIdx);
            choiceId_trl(isnan(choiceId_trl)) = [];
            angAcc_ln_exp(iSub, iB, iExp)     = length(find(choiceId_trl == 1)) / length(choiceId_trl);
        end
        % ------ different distractor numbers ------ 
        % angAcc_ln_dtr_exp = nan(subLen, BinL, 3, length(expList));
        for iB = 1 : BinL
            trlIdx       = (iB - 1) * nBin + 1 : iB * nBin;
            choiceId_trl = choiceId_ang(trlIdx);
            dtrNoCnt_trl = dtrNoCnt_iSub(trlIdx);
            dtrNoCnt_trl(isnan(choiceId_trl)) = [];
            choiceId_trl(isnan(choiceId_trl)) = [];
            for iDtr = 1 : 3
                angAcc_ln_dtr_exp(iSub, iB, iDtr, iExp) = length(find(choiceId_trl == 1 & dtrNoCnt_trl == iDtr)) / length(find(dtrNoCnt_trl == iDtr));
            end
        end

        % ------ within- vs. between-trans ------ 
        % angAcc_trans_ln_exp = nan(subLen, BinL, 2, length(expList));    % 2: within- vs. between transition
        for iB = 1 : BinL
            trlIdx         = (iB - 1) * nBin + 1 : iB * nBin;
            choiceId_trl   = choiceId_ang(trlIdx);
            transStyle_trl = transStyle(trlIdx);
            transStyle_trl(isnan(choiceId_trl)) = [];
            choiceId_trl(isnan(choiceId_trl)) = [];
            % ------ within-trans ------
            angAcc_trans_ln_exp(iSub, iB, 1, iExp) = length(find(choiceId_trl == 1 & transStyle_trl == 1)) / length(find(transStyle_trl == 1)); 
            % ------ between-trans------
            angAcc_trans_ln_exp(iSub, iB, 2, iExp) = length(find(choiceId_trl == 1 & transStyle_trl == 0)) / length(find(transStyle_trl == 0)); 
        end

        % ------ within- vs. between-trans for Random and Hamiltonian Walk ------ 
        % angAcc_trans_walk_ln_exp = nan(subLen, BinL, 2, 2, length(expList)); % first 2: within- vs. between transition; 2nd 2: random vs. hamiltonian walk
        for iB = 1 : BinL
            trlIdx         = (iB - 1) * nBin + 1 : iB * nBin;
            choiceId_trl   = choiceId_ang(trlIdx);
            transStyle_trl = transStyle(trlIdx);
            rndHam_Col_trl = rndHam_Col(trlIdx);
            transStyle_trl(isnan(choiceId_trl)) = [];
            rndHam_Col_trl(isnan(choiceId_trl)) = [];
            choiceId_trl(isnan(choiceId_trl))   = [];
            angAcc_trans_walk_ln_exp(iSub, iB, 1, 1, iExp) = length(find(choiceId_trl == 1 & transStyle_trl == 1 & rndHam_Col_trl == 1)) / length(find(transStyle_trl == 1 & rndHam_Col_trl == 1)); %% within & Random-walk
            angAcc_trans_walk_ln_exp(iSub, iB, 2, 1, iExp) = length(find(choiceId_trl == 1 & transStyle_trl == 0 & rndHam_Col_trl == 1)) / length(find(transStyle_trl == 0 & rndHam_Col_trl == 1)); %% between & Random-walk
            angAcc_trans_walk_ln_exp(iSub, iB, 1, 2, iExp) = length(find(choiceId_trl == 1 & transStyle_trl == 1 & rndHam_Col_trl == 2)) / length(find(transStyle_trl == 1 & rndHam_Col_trl == 2)); %% within & Hamiltonian-walk
            angAcc_trans_walk_ln_exp(iSub, iB, 2, 2, iExp) = length(find(choiceId_trl == 1 & transStyle_trl == 0 & rndHam_Col_trl == 2)) / length(find(transStyle_trl == 0 & rndHam_Col_trl == 2)); %% between & Hamiltonian-walk
        end

        %% accuracy of within- and between-transitions for Random and Hamiltonian Walk trials
        transStyle(isnan(choiceId_ang))   = [];
        rndHam_Col(isnan(choiceId_ang))   = [];
        nodesLabel(isnan(choiceId_ang))   = [];
        lureIn(isnan(choiceId_ang))       = [];
        choiceId_ang(isnan(choiceId_ang)) = [];
        acc_trans_walk(iSub, 1, 1, iExp) = length(find(choiceId_ang == 1 & transStyle == 1 & rndHam_Col == 1)) / length(find(transStyle == 1 & rndHam_Col == 1)); %% within & Random-walk
        acc_trans_walk(iSub, 2, 1, iExp) = length(find(choiceId_ang == 1 & transStyle == 0 & rndHam_Col == 1)) / length(find(transStyle == 0 & rndHam_Col == 1)); %% between & Random-walk
        acc_trans_walk(iSub, 1, 2, iExp) = length(find(choiceId_ang == 1 & transStyle == 1 & rndHam_Col == 2)) / length(find(transStyle == 1 & rndHam_Col == 2)); %% within & Hamiltonian-walk
        acc_trans_walk(iSub, 2, 2, iExp) = length(find(choiceId_ang == 1 & transStyle == 0 & rndHam_Col == 2)) / length(find(transStyle == 0 & rndHam_Col == 2)); %% between & Hamiltonian-walk

        %% calculate the trial proportions of within-cluster and between-clister transitions under Random and Hamiltonian Walk
        trlPro_trans_walk(iSub, 1, 1, iExp) = length(find(transStyle == 1 & rndHam_Col == 1)) / length(choiceId_ang);
        trlPro_trans_walk(iSub, 2, 1, iExp) = length(find(transStyle == 0 & rndHam_Col == 1)) / length(choiceId_ang);
        trlPro_trans_walk(iSub, 1, 2, iExp) = length(find(transStyle == 1 & rndHam_Col == 2)) / length(choiceId_ang);
        trlPro_trans_walk(iSub, 2, 2, iExp) = length(find(transStyle == 0 & rndHam_Col == 2)) / length(choiceId_ang);

        %% accuracy when lure stimulus exists
        % acc_trans_lure = zeros(subLen, 2, 2, length(expList)); % the 1st and 2nd '2' denote 'lure distractor exists vs. none' and 'transition from boundary node to within node vs. from boundary to boundary'
        % ------ lure stimulus exists in boundary-to-within trans ------ 
        acc_trans_lure(iSub, 1, 1, iExp) = length(find(choiceId_ang == 1 & nodesLabel == 1 & lureIn == 1 & transStyle == 1)) / length(find(nodesLabel == 1 & lureIn == 1 & transStyle == 1));
        % ------ lure stimulus exists in boundary-to-boundary trans ------ 
        acc_trans_lure(iSub, 1, 2, iExp) = length(find(choiceId_ang == 1 & nodesLabel == 1 & lureIn == 1 & transStyle ~= 1)) / length(find(nodesLabel == 1 & lureIn == 1 & transStyle ~= 1));
        % ------ no lure stimulus in boundary-to-within trans ------ 
        acc_trans_lure(iSub, 2, 1, iExp) = length(find(choiceId_ang == 1 & nodesLabel == 1 & lureIn == 0 & transStyle == 1)) / length(find(nodesLabel == 1 & lureIn == 0 & transStyle == 1));
        % ------ no lure stimulus in boundary-to-boundary trans ------ 
        acc_trans_lure(iSub, 2, 2, iExp) = length(find(choiceId_ang == 1 & nodesLabel == 1 & lureIn == 0 & transStyle ~= 1)) / length(find(nodesLabel == 1 & lureIn == 0 & transStyle ~= 1));

        % first 2:  with or without lure stimulus
        % second 2: boundary-to-within transition or bounary-to-boundary transition
        % third 3:  random or hamiltonian walk
        % acc_trans_walk_lure = zeros(subLen, 2, 2, 2, length(expList));
        % ------ lure stimulus exists ------
        % boundary-to-within trans, Random walk
        acc_trans_walk_lure(iSub, 1, 1, 1, iExp) = length(find(choiceId_ang == 1 & nodesLabel == 1 & lureIn == 1 & transStyle == 1 & rndHam_Col == 1)) / length(find(nodesLabel == 1 & lureIn == 1 & transStyle == 1 & rndHam_Col == 1));
        % boundary-to-boundary trans, Random walk
        acc_trans_walk_lure(iSub, 1, 2, 1, iExp) = length(find(choiceId_ang == 1 & nodesLabel == 1 & lureIn == 1 & transStyle ~= 1 & rndHam_Col == 1)) / length(find(nodesLabel == 1 & lureIn == 1 & transStyle ~= 1 & rndHam_Col == 1));
        % boundary-to-within trans, Hamiltonian walk
        acc_trans_walk_lure(iSub, 1, 1, 2, iExp) = length(find(choiceId_ang == 1 & nodesLabel == 1 & lureIn == 1 & transStyle == 1 & rndHam_Col == 2)) / length(find(nodesLabel == 1 & lureIn == 1 & transStyle == 1 & rndHam_Col == 2));
        % boundary-to-boundary trans, Hamiltonian walk
        acc_trans_walk_lure(iSub, 1, 2, 2, iExp) = length(find(choiceId_ang == 1 & nodesLabel == 1 & lureIn == 1 & transStyle ~= 1 & rndHam_Col == 2)) / length(find(nodesLabel == 1 & lureIn == 1 & transStyle ~= 1 & rndHam_Col == 2));
        % ------ no lure stimulus ------
        acc_trans_walk_lure(iSub, 2, 1, 1, iExp) = length(find(choiceId_ang == 1 & nodesLabel == 1 & lureIn == 0 & transStyle == 1 & rndHam_Col == 1)) / length(find(nodesLabel == 1 & lureIn == 0 & transStyle == 1 & rndHam_Col == 1));
        acc_trans_walk_lure(iSub, 2, 2, 1, iExp) = length(find(choiceId_ang == 1 & nodesLabel == 1 & lureIn == 0 & transStyle ~= 1 & rndHam_Col == 1)) / length(find(nodesLabel == 1 & lureIn == 0 & transStyle ~= 1 & rndHam_Col == 1));
        acc_trans_walk_lure(iSub, 2, 1, 2, iExp) = length(find(choiceId_ang == 1 & nodesLabel == 1 & lureIn == 0 & transStyle == 1 & rndHam_Col == 2)) / length(find(nodesLabel == 1 & lureIn == 0 & transStyle == 1 & rndHam_Col == 2));
        acc_trans_walk_lure(iSub, 2, 2, 2, iExp) = length(find(choiceId_ang == 1 & nodesLabel == 1 & lureIn == 0 & transStyle ~= 1 & rndHam_Col == 2)) / length(find(nodesLabel == 1 & lureIn == 0 & transStyle ~= 1 & rndHam_Col == 2));

        %% save data for the subsequent GLMM analysis
        % added by rxj @ 02/22/2024
        % data saved: choice (1-correct, 0-incorrect), trialNo, walkTypes,
        % transitionTypes
        

    end
    %%% merge all trials
    angAcc_exp(1 : subLen, :, iExp)             = angAcc;
    lenData_exp(1 : subLen, :, iExp)            = lenData;
    angAcc_dtr_exp(1 : subLen, :, :, iExp)      = angAcc_dtr;
    %%% only for the Random trials
    angAcc_Rand_exp(1 : subLen, :, iExp)        = angAcc_Rand;
    lenData_Rand_exp(1 : subLen, :, iExp)       = lenData_Rand;
    angAcc_dtr_Rand_exp(1 : subLen, :, :, iExp) = angAcc_dtr_Rand;
    %%% only for the Hamiltonian trials
    angAcc_Hami_exp(1 : subLen, :, iExp)        = angAcc_Hami;
    lenData_Hami_exp(1 : subLen, :, iExp)       = lenData_Hami;
    angAcc_dtr_Hami_exp(1 : subLen, :, :, iExp) = angAcc_dtr_Hami;

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

%% -------------SI figure: time curve by merging all trials with different distractors-------------
%% statistical testing in each time point
stat_mat = zeros(length(circle_list), 2, length(expList)); % 3 means three experiments
for iExp = 1 : length(expList)
    for iTime = 1 : length(circle_list)
        angAcc_exp_i = angAcc_exp(:, iTime, iExp);
        [h, p, ci, stats] = ttest(angAcc_exp_i, (1/2+1/3+1/4)/3, 'Tail', 'both');
        
        stat_mat(iTime, 1, iExp) = p;
        stat_mat(iTime, 2, iExp) = stats.tstat;
    end
end

%% SI figure: plotting the accuracy time curves
% angAcc_exp = zeros(subLen, length(circle_list), nExp); %% accuracy in each time point
LineSty  = '-';
colorTmp = colorSet(7, :);
for iExp = 1 : length(expList)
    figure('Position', [100 100 280 160]), clf;
    
    angAcc_exp_i = angAcc_exp(:, :, iExp);
    for iSub = 1 : size(angAcc_exp_i, 1)
        plot(circle_list, angAcc_exp_i(iSub, :), 'Color', [0.6, 0.6, 0.6], 'LineStyle', '-', 'LineWidth', 0.5); hold on;
    end
    [acc_avg, acc_sem] = Mean_and_Se(angAcc_exp_i, 1);
    shadedErrorBar(circle_list, acc_avg, acc_sem, {'Color', colorTmp, 'MarkerFaceColor', colorTmp, 'LineStyle', LineSty, 'LineWidth', 3}, 0.5); hold on;
    ylim([0, 1]);
    ylimit = ylim;
    xLoc = 0.99;
    % mark the significance
    [~, ~, ~, adj_p] = fdr_bh(squeeze(stat_mat(:, 1, iExp)), 0.05, 'pdep'); % method: 'dep', 'pdep'
    for iCir = 1 : length(circle_list)
        pval_j = adj_p(iCir);
        if pval_j < 0.05
            plot(circle_list(iCir), xLoc * ylimit(end), 'Marker', '.', 'MarkerSize', 4, 'Color', colorTmp, 'MarkerFaceColor', colorTmp, 'LineStyle', 'none'); hold on;
        end
    end
    plot(xlim, [(1/2+1/3+1/4)/3, (1/2+1/3+1/4)/3], 'k--', 'LineWidth', 1); hold on;
    if iExp == 1 || iExp == 3
        plot([0.8, 0.8], ylim, 'k--', 'LineWidth', 1); hold on;
    end
    ylim([0, 1]);
    axis xy;
    set(gca, 'LineWidth', 2);
    set(gca, 'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Arial');
    set(gca, 'XTick', 0 : 0.4 : 1.5, 'XTickLabel', '');
    set(gca, 'YTick', 0 : 0.5 : 1, 'YTickLabel', '');
    box off;
    ax = gca;
    %exportgraphics(ax, "Fig2e. learningCurve-LeafToHub.eps", "Resolution", 1000);
end

%% SI figure: choice accuracy trajectory for within- and between-transitions 
% %%% choice accuracy trajectory for within- and between-transitions
% angAcc_traj_trans_exp = nan(subLen, 2, length(circle_list), length(expList)); % 2: within and between-transitions
% %%% choice accuracy trajectory for within- and between-transitions with
% %%% final correct or incorrect response separation
% angAcc_traj_trans_cor_exp = nan(subLen, 4, length(circle_list), length(expList)); % Col1-2: final response correct vs. incorrect for the within-transitions; Col2: ~ for between-transitions
% dtrNo_traj_trans_cor_exp  = nan(subLen, 4, length(circle_list), length(expList)); % Quantifying the distractor number for each condition
dataFlg = 2;
if dataFlg == 1
    angAcc_traj_i = angAcc_traj_trans_exp;
elseif dataFlg == 2
    angAcc_traj_i = angAcc_traj_trans_cor_exp;
    dtrNo_traj_i  = dtrNo_traj_trans_cor_exp;
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
    markSize   = 4.5;
end
chance_i = mean([1/2, 1/3, 1/4]);
if dataFlg == 2
    chance_i_time = 1 ./ squeeze(nanmean(dtrNo_traj_i, 1)); % 4 * length(circle_list) * length(expList)
end
for iExp = 1 : 3
    disp(['---------- ', expList{iExp}, ' ----------']);
    figure('Position', [100 100 260 120]), clf;

    % ------ plot the chance level first ------
    if dataFlg == 2
        for ii = 1 : 4
            plot(circle_list, chance_i_time(ii, :, iExp), 'k--', 'LineWidth', 0.8); hold on;
        end
    end
    for iTb = 1 : 2     % within vs. between transitions
        if iTb == 1
            transWord = 'within';
            color_iTb = redGrad(2, :);
        elseif iTb == 2
            transWord = 'between';
            color_iTb = blueGrad(2, :);
        end
        if size(angAcc_traj_i, 2) == 2
            transData_traj_iTb = squeeze(angAcc_traj_i(:, iTb, :, iExp)); % subLen * length(circle_list) 
        elseif size(angAcc_traj_i, 2) == 4
            transData_traj_iTb = squeeze(angAcc_traj_i(:, (iTb - 1) * 2 + 1 : iTb * 2, :, iExp)); % subLen * 2* length(circle_list)
        end
        [accAvg_rOh, accSem_rOh] = Mean_and_Se(transData_traj_iTb, 1);  
        accAvg_rOh = squeeze(accAvg_rOh);
        accSem_rOh = squeeze(accSem_rOh);

        lineLen = size(accAvg_rOh, 1);
        for ilr = 1 : lineLen % final correct vs. incorrect response
            if ilr == 1
                LineSty = '-';
            elseif ilr == 2
                LineSty = ':';
            end
            shadedErrorBar(circle_list, accAvg_rOh(ilr, :), accSem_rOh(ilr, :), {'Color', color_iTb, 'MarkerFaceColor', color_iTb, 'LineStyle', LineSty, 'LineWidth', 2}, 0.5); hold on;
        end
    end
    xlim([0, 0.8]);
    ylim([0, 1]);
    if dataFlg == 1
        plot(xlim, [chance_i, chance_i], 'k--', 'LineWidth', 0.8); hold on;
    end
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
        set(gca, 'XTick', 0 : 0.4 : 0.8, 'XTickLabel', '');
        set(gca, 'YTick', 0 : 0.5 : 1, 'YTickLabel', '');
    end
    box off;
end

%% SI figure: learning curve for different distractors across trial bins (100 trials/bin)
% angAcc_ln_dtr_exp = nan(subLen, BinL, 3, length(expList));
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
% ------ statistical tests for each bin ------
statsMat_dtr_exp = nan(length(dtrNums), BinL, 2, length(expList)); % 2: p and tstats
for iExp = 1 : length(expList)
    for iDt = 1 : length(dtrNums)
        ref_iDt = 1 / (iDt + 1);
        for iB = 1 : BinL
            angAcc_iB = squeeze(angAcc_ln_dtr_exp(:, iB, iDt, iExp));
            % one-sample test
            [h, p, ci, stats] = ttest(angAcc_iB, ref_iDt, 'Tail', 'right');
            statsMat_dtr_exp(iDt, iB, 1, iExp) = p;
            statsMat_dtr_exp(iDt, iB, 2, iExp) = stats.tstat;
        end
    end
end
LineStys = {'-', '-.', ':'};
for iExp = 1 : length(expList)
    figure('Position', [100 100 260 160]), clf;
    hold on;
    xlim([0.5, BinL+0.5]);
    ylim([0.2, 0.8]);
    plot(xlim, [1/2, 1/2], 'Color', [0.6, 0.6, 0.6], 'LineStyle', LineStys{1}, 'LineWidth', 0.8); hold on;
    plot(xlim, [1/3, 1/3], 'Color', [0.6, 0.6, 0.6], 'LineStyle', LineStys{2}, 'LineWidth', 0.8); hold on;
    plot(xlim, [1/4, 1/4], 'Color', [0.6, 0.6, 0.6], 'LineStyle', LineStys{3}, 'LineWidth', 0.8); hold on;
    for iDt = 1 : length(dtrNums) 
        angAcc_iDt = squeeze(angAcc_ln_dtr_exp(:, :, iDt, iExp));
        [accAvg, accSem] = Mean_and_Se(angAcc_iDt, 1);
        errorbar(1 : 1 : BinL, accAvg, accSem, 'Color', blueGrad(iDt, :), 'LineStyle', LineStys{iDt}, 'LineWidth', errLineWid); hold on;

        [~, ~, ~, adj_p] = fdr_bh(squeeze(statsMat_dtr_exp(iDt, :, 1, iExp)), 0.05, 'pdep'); % method: 'dep', 'pdep'
        for iB = 1 : BinL
            if adj_p(iB) < 0.05
                markFace = blueGrad(iDt, :);
            else
                markFace = [1, 1, 1];
            end
            plot(iB, accAvg(iB), 'Marker', 'o', 'MarkerSize', markSize, 'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', markFace, 'LineStyle', '-'); hold on;
        end
    end
    xlim([0.5, BinL+0.5]);
    ylim([0.2, 0.8]);
    if figKey == 0
        % ------For presentation------
        set(gca, 'LineWidth', 2);
        set(gca, 'FontSize', 15, 'FontWeight', 'bold', 'FontName', 'Arial');
        set(gca, 'XTick', 1 : 1 : BinL, 'XTickLabel', 1 : 1 : BinL);
        set(gca, 'YTick', 0.2 : 0.2 : 0.8, 'YTickLabel', 0.2 : 0.2 : 0.8);
    elseif figKey == 1
        % ------For Adobe Illustrator------
        set(gca, 'LineWidth', 0.8);
        set(gca, 'FontSize', 10, 'FontWeight', 'bold', 'FontName', 'Arial');
        set(gca, 'XTick', 1 : 1 : BinL, 'XTickLabel', '');
        set(gca, 'YTick', 0.2 : 0.2 : 0.8, 'YTickLabel', '');
    end
    box off;
end

%% BhevaiorPaper Figure xx: learning curve for within and between transitions across trial bins (100 trials/bin)
% angAcc_trans_ln_exp = nan(subLen, BinL, 2, length(expList));    % 2: within- vs. between transition
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
nWalks = 2; % Random and Hamiltonian Walk
nTrans = 2; % within and Between-cluster transition
% ------ statistical tests for each bin ------
statsMat_trans_exp = nan(nTrans, BinL, 2, length(expList)); % 2: p and tstats
ref_trans = (1/2+1/3+1/4)/3;
for iExp = 1 : length(expList)
    for iTs = 1 : nTrans % within- vs. between-trans
        for iB = 1 : BinL
            angAcc_iB = squeeze(angAcc_trans_ln_exp(:, iB, iTs, iExp));
            % one-sample test
            [h, p, ci, stats] = ttest(angAcc_iB, ref_trans, 'Tail', 'right');
            statsMat_trans_exp(iTs, iB, 1, iExp) = p;
            statsMat_trans_exp(iTs, iB, 2, iExp) = stats.tstat;
        end
    end
end
for iExp = 1 : length(expList)
    figure('Position', [100 100 260 160]), clf;
    hold on;
    xlim([0.5, BinL+0.5]);
    ylim([0.2, 0.8]);
    plot(xlim, [(1/2+1/3+1/4)/3, (1/2+1/3+1/4)/3], 'Color', [0.6, 0.6, 0.6], 'LineStyle', LineStys{1}, 'LineWidth', 0.8); hold on;

    for iTs = 1 : nTrans % within- vs. between-trans
        angAcc_ii = squeeze(angAcc_trans_ln_exp(:, :, iTs, iExp));
        [accAvg, accSem] = Mean_and_Se(angAcc_ii, 1);
        colorTmp = colorSet(iTs, :);
        errorbar(1 : 1 : BinL, accAvg, accSem, 'Color', colorTmp, 'LineStyle', '-', 'LineWidth', errLineWid); hold on;

        [~, ~, ~, adj_p] = fdr_bh(squeeze(statsMat_trans_exp(iTs, :, 1, iExp)), 0.05, 'pdep'); % method: 'dep', 'pdep'
        for iB = 1 : BinL
            if adj_p(iB) < 0.05
                markFace = colorTmp;
            else
                markFace = [1, 1, 1];
            end
            plot(iB, accAvg(iB), 'Marker', 'o', 'MarkerSize', markSize, 'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', markFace, 'LineStyle', '-'); hold on;
        end
    end
    xlim([0.5, BinL+0.5]);
    ylim([0.2, 0.8]);
    if figKey == 0
        % ------For presentation------
        set(gca, 'LineWidth', 2);
        set(gca, 'FontSize', 15, 'FontWeight', 'bold', 'FontName', 'Arial');
        set(gca, 'XTick', 1 : 1 : BinL, 'XTickLabel', 1 : 1 : BinL);
        set(gca, 'YTick', 0.2 : 0.2 : 0.8, 'YTickLabel', 0.2 : 0.2 : 0.8);
    elseif figKey == 1
        % ------For Adobe Illustrator------
        set(gca, 'LineWidth', 0.8);
        set(gca, 'FontSize', 10, 'FontWeight', 'bold', 'FontName', 'Arial');
        set(gca, 'XTick', 1 : 1 : BinL, 'XTickLabel', '');
        set(gca, 'YTick', 0.2 : 0.2 : 0.8, 'YTickLabel', '');
    end
    box off;
end

%% BehaviorPaper Figure xx: learning curve for within and between transitions in Random and Hamiltonian Walk across trial bins (100 trials/bin)
% angAcc_trans_walk_ln_exp = nan(subLen, BinL, 2, 2, length(expList)); % first 2: within- vs. between transition; 2nd 2: random vs. hamiltonian walk
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
nWalks = 2; % Random and Hamiltonian Walk
nTrans = 2; % within and Between-cluster transition
for iExp = 1 : length(expList)
    figure('Position', [100 100 260 160]), clf;
    hold on;
    xlim([0.5, BinL+0.5]);
    ylim([0.2, 0.8]);
    plot(xlim, [(1/2+1/3+1/4)/3, (1/2+1/3+1/4)/3], 'Color', [0.6, 0.6, 0.6], 'LineStyle', LineStys{1}, 'LineWidth', 0.8); hold on;

    for i_rOh = 1 : nWalks   % random walk vs. hamiltonian walk
        for iTs = 1 : nTrans % within- vs. between-trans
            angAcc_ii = squeeze(angAcc_trans_walk_ln_exp(:, :, iTs, i_rOh, iExp));
            [accAvg, accSem] = Mean_and_Se(angAcc_ii, 1);
            if i_rOh == 1
                colorTmp = colorSet(iTs, :);
            elseif i_rOh == 2
                colorTmp = 0.5 * colorSet(iTs, :) + 0.5 * [1, 1, 1];
            end
            errorbar(1 : 1 : BinL, accAvg, accSem, 'Color', colorTmp, 'LineStyle', '-', 'LineWidth', errLineWid); hold on;
            for iB = 1 : BinL
                plot(iB, accAvg(iB), 'Marker', 'o', 'MarkerSize', markSize, 'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', colorTmp, 'LineStyle', '-'); hold on;
            end
        end
    end
    xlim([0.5, BinL+0.5]);
    ylim([0.2, 0.8]);
    if figKey == 0
        % ------For presentation------
        set(gca, 'LineWidth', 2);
        set(gca, 'FontSize', 15, 'FontWeight', 'bold', 'FontName', 'Arial');
        set(gca, 'XTick', 1 : 1 : BinL, 'XTickLabel', 1 : 1 : BinL);
        set(gca, 'YTick', 0.2 : 0.2 : 0.8, 'YTickLabel', 0.2 : 0.2 : 0.8);
    elseif figKey == 1
        % ------For Adobe Illustrator------
        set(gca, 'LineWidth', 0.8);
        set(gca, 'FontSize', 10, 'FontWeight', 'bold', 'FontName', 'Arial');
        set(gca, 'XTick', 1 : 1 : BinL, 'XTickLabel', '');
        set(gca, 'YTick', 0.2 : 0.2 : 0.8, 'YTickLabel', '');
    end
    box off;
end

%% BehaviorPaper Figure xx: accuracy for the within- and between-transitions in Random and Hamiltonian Walk trials separately
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

for iExp = 1 : 3
    disp(['---------- ', expList{iExp}, ' ----------']);
    figure('Position', [100 100 260 120]), clf;

    acc_trans_walk_iExp = acc_trans_walk(:, :, :, iExp);
    [accAvg_rOh, accSem_rOh] = Mean_and_Se(acc_trans_walk_iExp, 1);
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
        plot(barPos_i, acc_trans_walk_iExp(:, :, i_rOh), 'Color', [0.6, 0.6, 0.6], 'LineStyle', '-', 'LineWidth', indvLineW); hold on;
        plot(barPos_i, accAvg_rOh(:, i_rOh), 'Color', [0, 0, 0], 'LineStyle', '-', 'LineWidth', errLineWid); hold on;
        for i_tp = 1 : 2 % within vs. between-transition
            if i_rOh == 1
                colorTmp = colorSet(i_tp, :);
            elseif i_rOh == 2
                colorTmp = 0.5 * colorSet(i_tp, :) + 0.5 * [1, 1, 1];
            end
            errorbar(barPos_i(i_tp), accAvg_rOh(i_tp, i_rOh), accSem_rOh(i_tp, i_rOh), 'Color', 'k', 'LineStyle', 'none', 'LineWidth', errLineWid); hold on;
            plot(barPos_i(i_tp), accAvg_rOh(i_tp, i_rOh), 'Marker', 'o', 'MarkerSize', markSize, 'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', colorTmp, 'LineStyle', '-'); hold on;
        end
        % ------ Statistical tests ------
        disp(['======== ', walkWord, '-within vs. between trans ========']);
        [h, p, ci, stats] = ttest(acc_trans_walk_iExp(:, 1, i_rOh), acc_trans_walk_iExp(:, 2, i_rOh))
    end
    xlim([0.6, 2.6]);
    ylim([0.05, 0.8]);
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
    box off;

end

%% paired-samplle t-test for within- and between-transitions
walkNo = 2;
pairT_p = nan(length(expList), walkNo); % walkNo: random and hamiltonian walk
pairT_t = nan(length(expList), walkNo);
for iExp = 1 : length(expList)
    for i_rOh = 1 : 2 % Random vs. Hamiltonian Walk
        acc_trans_walk_i  = acc_trans_walk(:, :, i_rOh, iExp);
        [h, p, ci, stats] = ttest(acc_trans_walk_i(:, 1), acc_trans_walk_i(:, 2));
        pairT_p(iExp, i_rOh) = p;
        pairT_t(iExp, i_rOh) = stats.tstat;
    end
end

%% one-sample ttest: comparison with chance level
transNo = 2; % within- and between-trans
walkNo  = 2; % random and hamiltonian walk
oneSampT_stat = nan(walkNo, transNo, 2, length(expList)); % 2: p and t-stats
chance_i = mean([1/2, 1/3, 1/4]);
for iExp = 1 : length(expList)
    for i_rOh = 1 : 2    % Random vs. Hamiltonian Walk
        for i_tp = 1 : 2 % Within vs. Between transition
            acc_trans_walk_i  = acc_trans_walk(:, i_tp, i_rOh, iExp);
            [h, p, ci, stats] = ttest(acc_trans_walk_i, chance_i);
            oneSampT_stat(i_rOh, i_tp, 1, iExp) = p;
            oneSampT_stat(i_rOh, i_tp, 2, iExp) = stats.tstat;
        end
    end
end

%% with lure distractors vs. without: one-hot metric
% acc_trans_lure      = zeros(subLen, 2, 2, length(expList)); % the 1st and 2nd '2' denote 'lure distractor exists vs. none' and 'transition from boundary node to within node vs. from boundary to boundary'
% acc_trans_walk_lure = zeros(subLen, 2, 2, 2, length(expList));
% first 2:  with or without lure stimulus
% second 2: boundary-to-within transition or bounary-to-boundary transition
% third 3:  random or hamiltonian walk

flg = 3;
if flg == 1
    ACCmat = acc_trans_lure; %% mixture of random and hamiltonian path
elseif flg == 2
    ACCmat = squeeze(acc_trans_walk_lure(:, :, :, 1, :)); %% only random path
elseif flg == 3
    ACCmat = squeeze(acc_trans_walk_lure(:, :, :, 2, :)); %% only hamiltonian path
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
chance_i = mean([1/2, 1/3, 1/4]);
for iExp = 1 : 3
    disp(['---------- ', expList{iExp}, ' ----------']);
    figure('Position', [100 100 260 120]), clf;

    ACCmat_iExp = ACCmat(:, :, :, iExp);
    [accAvg_rOh, accSem_rOh] = Mean_and_Se(ACCmat_iExp, 1);
    accAvg_rOh = squeeze(accAvg_rOh);
    accSem_rOh = squeeze(accSem_rOh);
    barPos = [1, 1.5; 1.7, 2.2];
    for iTb = 1 : 2     % 'transition from boundary node to within node vs. from boundary to boundary'
        if iTb == 1
            transWord = 'boundary-to-within';
        elseif iTb == 2
            transWord = 'boundary-to-boundary';
        end
        barPos_i = barPos(iTb, :);
        %%% line plot
        plot(barPos_i, ACCmat_iExp(:, :, iTb), 'Color', [0.6, 0.6, 0.6], 'LineStyle', '-', 'LineWidth', indvLineW); hold on;
        plot(barPos_i, accAvg_rOh(:, iTb), 'Color', [0, 0, 0], 'LineStyle', '-', 'LineWidth', errLineWid); hold on;
        for ilr = 1 : 2 % 'lure distractor exists vs. none'
            if ilr == 1
                colorTmp = [0, 0, 0]; % with lure distractor
            elseif ilr == 2
                colorTmp = [1, 1, 1]; % without lure distractor
            end
            errorbar(barPos_i(ilr), accAvg_rOh(ilr, iTb), accSem_rOh(ilr, iTb), 'Color', 'k', 'LineStyle', 'none', 'LineWidth', errLineWid); hold on;
            plot(barPos_i(ilr), accAvg_rOh(ilr, iTb), 'Marker', 'o', 'MarkerSize', markSize, 'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', colorTmp, 'LineStyle', '-'); hold on;
        end
        % ------ Statistical tests ------
        disp(['======== ', transWord, ': with vs. without lure stimulus ========']);
        [h, p, ci, stats] = ttest(ACCmat_iExp(:, 1, iTb), ACCmat_iExp(:, 2, iTb))
    end
    xlim([0.6, 2.6]);
    ylim([0.05, 0.8]);
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
    box off;
end

%% with lure distractors vs. without: trajectory
% acc_trans_traj_lure      = zeros(subLen, length(circle_list), 2, 2, length(expList)); % the 1st and 2nd '2' denote 'lure distractor exists vs. none' and 'transition from boundary node to within node vs. from boundary to boundary'
% acc_trans_walk_traj_lure = zeros(subLen, length(circle_list), 2, 2, 2, length(expList));
% first 2:  with or without lure stimulus
% second 2: boundary-to-within transition or bounary-to-boundary transition
% third 3:  random or hamiltonian walk

flg = 1;
if flg == 1
    ACCmat_traj = acc_trans_traj_lure; %% mixture of random and hamiltonian path
elseif flg == 2
    ACCmat_traj = squeeze(acc_trans_walk_traj_lure(:, :, :, :, 1, :)); %% only random path
elseif flg == 3
    ACCmat_traj = squeeze(acc_trans_walk_traj_lure(:, :, :, :, 2, :)); %% only hamiltonian path
end

%% BehavioralPaper, Figure xx: influence of lure stimulus on choice accuracy trajectory
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
chance_i = mean([1/2, 1/3, 1/4]);

for iExp = 1 : 3
    disp(['---------- ', expList{iExp}, ' ----------']);
    figure('Position', [100 100 260 120]), clf;

    ACCmat_iExp = ACCmat_traj(:, :, :, :, iExp);
    [accAvg_rOh, accSem_rOh] = Mean_and_Se(ACCmat_iExp, 1);
    accAvg_rOh = squeeze(accAvg_rOh); % length(circle_list) * 2 * 2
    accSem_rOh = squeeze(accSem_rOh);
    for iTb = 1 : 2     % 'transition from boundary node to within node vs. from boundary to boundary'
        if iTb == 1
            transWord = 'boundary-to-within';
            color_iTb = redGrad(2, :);
        elseif iTb == 2
            transWord = 'boundary-to-boundary';
            color_iTb = blueGrad(2, :);
        end
        for ilr = 1 : 2 % 'lure distractor exists vs. none'
            if ilr == 1
                LineSty = '-';
            elseif ilr == 2
                LineSty = ':';
            end
            shadedErrorBar(circle_list, accAvg_rOh(:, ilr, iTb), accSem_rOh(:, ilr, iTb), {'Color', color_iTb, 'MarkerFaceColor', color_iTb, 'LineStyle', LineSty, 'LineWidth', 2}, 0.5); hold on;
        end
    end
    xlim([0, 0.8]);
    ylim([0, 1]);
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
        set(gca, 'YTick', 0 : 0.5 : 1, 'YTickLabel', '');
    end
    box off;
end


