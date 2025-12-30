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

%% =================== Model fitting procedure ===================
%% Experiment 1:
% I didn't regress out the confounding variables from RTs.
circle_list = 0 : 1/60 : 1.5;
states      = 15;
nFit        = 100;
expList     = {'ImplicitExp', 'ExplicitExp', 'ImplicitRandExp'}; % {'ImplicitExp', 'ExplicitExp', 'ImplicitRandExp', 'ImplicitClusterExp', 'ImplicitRandClusterExp'};
folder      = '/Users/renxiangjuan/Nextcloud/MATLAB/MyExperiment/HierarchicalCluster'; %'/Volumes/HierarchicalCluster';
%folder      = '/Users/cdlab/Nextcloud/MATLAB/MyExperiment/HierarchicalCluster'; %pwd;
%folder          ='/home/ubuntu/Xiangjuan/MyExperiment/HierarchicalCluster';

subLen = 24;
%%% merge all trials
angAcc_exp     = zeros(subLen, length(circle_list), length(expList)); %% accuracy in each time point
lenData_exp    = zeros(subLen, length(circle_list), length(expList)); %% data length in each time point
angAcc_dtr_exp = zeros(subLen, length(circle_list), 3, length(expList));
%%% only for the Random trials
angAcc_Rand_exp     = zeros(subLen, length(circle_list), length(expList)); %% accuracy in each time point
lenData_Rand_exp    = zeros(subLen, length(circle_list), length(expList)); %% data length in each time point
angAcc_dtr_Rand_exp = zeros(subLen, length(circle_list), 3, length(expList));
%%% only for the Hamiltonian trials
angAcc_Hami_exp     = zeros(subLen, length(circle_list), length(expList)); %% accuracy in each time point
lenData_Hami_exp    = zeros(subLen, length(circle_list), length(expList)); %% data length in each time point
angAcc_dtr_Hami_exp = zeros(subLen, length(circle_list), 3, length(expList));
%%% accuracy for the within and between transitions for Random and Hamiltonian Walk separately
acc_trans_walk = zeros(subLen, 2, 2, length(expList)); % the 1st and 2nd '2' denote 'within vs. between cluster transition' and 'random and hamiltonian' walk
%%% trial proportions of within-cluster and between-cluster transition under Random and Hamiltonian Walk
trlPro_trans_walk = zeros(subLen, 2, 2, length(expList)); 
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
        
    elseif isequal(ExpWord, 'ImplicitClusterExp')
        bhvDataDir = [folder, '/FormalExp-ClusterComponentLearning-Results/'];
        subj_list = {'gyl_1_f_21', 'yrl_2_f_20', 'xwh_3_m_19', 'lz_4_m_20', 'mrs_5_f_24', 'wy_6_m_24', 'ljx_7_f_26', 'hhn_8_f_19', ...
                     'rhj_9_f_27', 'lya_10_f_20', 'lfy_11_m_24'};
        subjLab   = {'gyl1', 'yrl2', 'xwh3', 'lz4', 'mrs5', 'wy6', 'ljx7', 'hhn8', ...
                     'rhj9', 'lya10', 'lfy11'};
        tInBlc = [200, 2*200, 3*200, 3*200+100, rndTrial + (200 : 200 : HamTrial)];
        
    elseif isequal(ExpWord, 'ImplicitRandClusterExp')
        bhvDataDir = [folder, '/FormalExp-ImplicitRandom-ClusterComponentLearning-Results/'];
        subj_list = {'xy_1_f_21', 'xly_2_f_21', 'wyt_3_f_23', 'zty_4_m_18', 'lyn_5_m_19', 'wyq_6_f_19', 'sy_7_f_22', 'wyh_8_f_27', 'qbj_9_f_25', 'lt_10_m_23', ...
                     'gjy_11_f_23', 'wyh_12_m_23', 'zwq_13_m_22'}; % , 'wc_14_f_23'
        subjLab   = {'xy1', 'xly2', 'wyt3', 'zty4', 'lyn5', 'wyq6', 'sy7', 'wyh8', 'qbj9', 'lt10', ...
                     'gjy11', 'wyh12', 'zwq13'}; % , 'wc14' 
        tInBlc = [200, 2*200, 3*200, 3*200+100, rndTrial + (200 : 200 : HamTrial)];
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
    
    %% fitting for single subject
    %%% merge all trials
    angAcc  = zeros(subLen, length(circle_list)); %% accuracy in each time point
    lenData = zeros(subLen, length(circle_list)); %% data length in each time point
    angAcc_dtr = zeros(subLen, length(circle_list), 3); % 3 means distractor number categories
    %%% only for the Random trials
    angAcc_Rand  = zeros(subLen, length(circle_list)); %% accuracy in each time point
    lenData_Rand = zeros(subLen, length(circle_list)); %% data length in each time point
    angAcc_dtr_Rand = zeros(subLen, length(circle_list), 3); % 3 means distractor number categories
    %%% only for the Hamiltonian trials
    angAcc_Hami  = zeros(subLen, length(circle_list)); %% accuracy in each time point
    lenData_Hami = zeros(subLen, length(circle_list)); %% data length in each time point
    angAcc_dtr_Hami = zeros(subLen, length(circle_list), 3); % 3 means distractor number categories
    for iSub = 1 : subLen
        disp([ExpWord, '-subj', num2str(iSub)]);
        %%
        subID  = subjLab{iSub};
        subjBv = subj_listBv{iSub};
        subjDir = [bhvDataDir, subjBv, '/'];
        if isequal(ExpWord, 'ImplicitClusterExp') || isequal(ExpWord, 'ImplicitRandClusterExp')
            load([subjDir, subjBv, 'clusterResult_reOrg.mat'], 'clusterResult');
        else
            load([subjDir, subjBv, 'clusterResult_Blocks.mat'], 'clusterResult');
        end
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
        
        %%% mouse trajectory trials
        load([subjDir, subjBv,  'respYes_trials_blc.mat'], 'respYes_trials_blc');
        load([subjDir, subjBv,  'mouseTraj_trials_blc.mat'], 'mouseTraj_trials_blc');
        if isequal(ExpWord, 'ImplicitRandClusterExp')
            load([subjDir, subjBv, 'trialsInBlc.mat'], 'trialsInBlc');
        end
        
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
            choiceId_i = choiceId(:, iTp);
            dtrNoCnt_iSub_i = dtrNoCnt_iSub;
            dtrNoCnt_iSub_i(isnan(choiceId_i)) = [];
            choiceId_i(isnan(choiceId_i))     = [];
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
            choiceId_i = choiceId(used_Trl, iTp);
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
            choiceId_i = choiceId(used_Trl, iTp);
            dtrNoCnt_iSub_i = dtrNoCnt_iSub(used_Trl);
            dtrNoCnt_iSub_i(isnan(choiceId_i)) = [];
            choiceId_i(isnan(choiceId_i))     = [];
            angAcc_Hami(iSub, iTp)  = length(find(choiceId_i == 1)) / length(choiceId_i);
            lenData_Hami(iSub, iTp) = length(choiceId_i);

            for iDtr = 1 : 3
                angAcc_dtr_Hami(iSub, iTp, iDtr) = length(find(choiceId_i == 1 & dtrNoCnt_iSub_i == iDtr)) / length(find(dtrNoCnt_iSub_i == iDtr));
            end
        end

        %% accuracy of within- and between-transitions for Random and Hamiltonian Walk trials
        transStyle(isnan(choiceId_ang))   = [];
        rndHam_Col(isnan(choiceId_ang))   = [];
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

        %% save data for the subsequent GLMM analysis
        % added by rxj @ 02/22/2024
        % data saved: choice (1-correct, 0-incorrect), trialNo, walkTypes,
        % transitionTypes
        

    end
    %%% merge all trials
    angAcc_exp(1 : subLen, :, iExp)  = angAcc;
    lenData_exp(1 : subLen, :, iExp) = lenData;
    angAcc_dtr_exp(1 : subLen, :, :, iExp) = angAcc_dtr;
    %%% only for the Random trials
    angAcc_Rand_exp(1 : subLen, :, iExp)  = angAcc_Rand;
    lenData_Rand_exp(1 : subLen, :, iExp) = lenData_Rand;
    angAcc_dtr_Rand_exp(1 : subLen, :, :, iExp) = angAcc_dtr_Rand;
    %%% only for the Hamiltonian trials
    angAcc_Hami_exp(1 : subLen, :, iExp)  = angAcc_Hami;
    lenData_Hami_exp(1 : subLen, :, iExp) = lenData_Hami;
    angAcc_dtr_Hami_exp(1 : subLen, :, :, iExp) = angAcc_dtr_Hami;

end


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
colorComp = colorSets([1, 3, 5], :);

%% -------------Part 1: merge all trials with different distractors-------------
%% statistical testing in each time point
stat_mat = zeros(length(circle_list), 2, length(expList)); % 3 means three experiments
for iExp = 1 : length(expList)
    for iTime = 1 : length(circle_list)
        angAcc_exp_i = angAcc_exp(:, iTime, iExp);
        [h, p, ci, stats] = ttest(angAcc_exp_i, (1/2+1/3+1/4)/3, 'Tail', 'right');
        
        stat_mat(iTime, 1, iExp) = p;
        stat_mat(iTime, 2, iExp) = stats.tstat;
    end
end

%% plotting the accuracy time curves
% angAcc  = zeros(subLen, length(circle_list));
figure('Position', [100 100 1000 300]), clf;
ha = tight_subplot(1, length(expList), [.05 .05], [.1 .02], [.05 .05]);
ha_i = 1;
for iExp = 1 : length(expList)
    angAcc_exp_i  = angAcc_exp(:, :, iExp);
    lenData_exp_i = lenData_exp(:, :, iExp);
    [acc_avg, acc_sem] = Mean_and_Se(angAcc_exp_i, 1);
    
    axes(ha(ha_i));
    %% plotting the curve for individual subject
    for iSub = 1 : subLen
        plot(circle_list, angAcc_exp_i(iSub, :), 'Color', [0.4, 0.4, 0.4], 'LineStyle', '-', 'LineWidth', 1); hold on;
    end
    shadedErrorBar(circle_list, acc_avg, acc_sem, {'Color', 'r', 'MarkerFaceColor', 'r', 'LineStyle', '-', 'LineWidth', 2}, 0.5); hold on;
    ylim([0, 1]);
    ylimit = ylim;
    xLoc = 0.98;
    %%
    [~, ~, ~, adj_p] = fdr_bh(squeeze(stat_mat(:, 1, iExp)), 0.05, 'pdep'); % method: 'dep', 'pdep'
    for iCir = 1 : length(circle_list)
        pval_j = adj_p(iCir);
        if pval_j < 0.01
            plot(circle_list(iCir), xLoc * ylimit(end), 'Marker', '.', 'MarkerSize', 8, 'Color', colorSets(iExp, :), 'MarkerFaceColor', colorSets(iExp, :), 'LineStyle', 'none', 'LineWidth', 1.5); hold on;
        elseif pval_j < 0.05 && pval_j >= 0.01
            plot(circle_list(iCir), xLoc * ylimit(end), 'Marker', '>', 'MarkerSize', 4, 'Color', colorSets(iExp, :), 'MarkerFaceColor', colorSets(iExp, :), 'LineStyle', 'none', 'LineWidth', 1.5); hold on;
        end
    end
    if iExp == 1 || iExp == 3
        plot([0.8, 0.8], ylim, 'k:', 'LineWidth', 1); hold on;
    end
    plot(xlim, [(1/2+1/3+1/4)/3, (1/2+1/3+1/4)/3], 'k--', 'LineWidth', 1); hold on;
    axis xy;
    set(gca, 'LineWidth', 2);
    set(gca, 'FontSize', 15, 'FontWeight', 'bold', 'FontName', 'Arial');
    box off;
    
    ha_i = ha_i + 1;
end
%% plotting the accuracy from 3 experiments together
% added by rxj @ 06/14/2022
% angAcc_exp = zeros(subLen, length(circle_list), length(expList)); %% accuracy in each time point
figure('Position', [100 100 600 220]), clf;
test_name  = 'oneSampT';
xLoc = 0.99;
for iExp = 1 : length(expList)
    angAcc_times = angAcc_exp(:, :, iExp);

    %% statistical test
    %%% one-sample t test
    tail_id = 'right';
    stat_mat = zeros(length(circle_list), 2); % 2: p value and t score
    for iTime = 1 : length(circle_list)
        angAcc_exp_i = angAcc_times(:, iTime);
        [h, p, ci, stats] = ttest(angAcc_exp_i, (1/2+1/3+1/4)/3, 'Tail', tail_id);
        stat_mat(iTime, 1) = p;
        stat_mat(iTime, 2) = stats.tstat;
    end

    %% plotting
    [beta_Avg, beta_Sem] = Mean_and_Se(angAcc_times, 1);
    shadedErrorBar(circle_list, beta_Avg, beta_Sem, {'Color', colorComp(iExp, :), 'MarkerFaceColor', colorComp(iExp, :), 'LineStyle', '-', 'LineWidth', 3}, 0.5); hold on;
    ylim([0, 1]);
    ylimit = ylim;
    %% FDR correction
    % stat_mat = zeros(length(circle_list), 2); % 2: p value and t score
    [~, ~, ~, adj_p] = fdr_bh(stat_mat(:, 1), 0.05, 'pdep'); % method: 'dep', 'pdep'
    for iTime = 1 : length(circle_list)
        p_i = adj_p(iTime);
        if p_i < 0.01
            plot(circle_list(iTime), xLoc * ylimit(end), 'Marker', '*', 'MarkerSize', 4, 'Color', colorComp(iExp, :)); hold on;
        elseif p_i < 0.05 && p_i >= 0.01
            plot(circle_list(iTime), xLoc * ylimit(end), 'Marker', '.', 'MarkerSize', 10, 'Color', colorComp(iExp, :), 'MarkerFaceColor', colorComp(iExp, :), 'LineStyle', 'none', 'LineWidth', 1.5); hold on;
        end
    end
    xLoc = xLoc - 0.05;
end
%xlim([0, 0.8]);
plot(xlim, [(1/2+1/3+1/4)/3, (1/2+1/3+1/4)/3], 'k--', 'LineWidth', 1); hold on;
plot([0.8, 0.8], ylim, 'k:', 'LineWidth', 1); hold on;
set(gca, 'LineWidth', 2);
set(gca, 'FontSize', 18, 'FontWeight', 'bold', 'FontName', 'Arial');
box off;

%% plotting the accuracy from 3 experiments together: separately for training and testing
% added by rxj @ 06/14/2022
%angAcc_Rand_exp = zeros(subLen, length(circle_list), length(expList)); %% accuracy in each time point
%angAcc_Hami_exp = zeros(subLen, length(circle_list), length(expList)); %% accuracy in each time point
for iPath = 1 : 2
    if iPath == 1
        angAcc_iPath = angAcc_Rand_exp;
    elseif iPath == 2
        angAcc_iPath = angAcc_Hami_exp;
    end
    figure('Position', [100 100 400 220]), clf;
    test_name  = 'oneSampT';
    xLoc = 0.99;
    for iExp = 1 : length(expList)
        angAcc_times = angAcc_iPath(:, :, iExp);
        
        %% statistical test
        %%% one-sample t test
        tail_id = 'right';
        stat_mat = zeros(length(circle_list), 2); % 2: p value and t score
        for iTime = 1 : length(circle_list)
            angAcc_exp_i = angAcc_times(:, iTime);
            [h, p, ci, stats] = ttest(angAcc_exp_i, (1/2+1/3+1/4)/3, 'Tail', tail_id);
            stat_mat(iTime, 1) = p;
            stat_mat(iTime, 2) = stats.tstat;
        end
        
        %% plotting
        [beta_Avg, beta_Sem] = Mean_and_Se(angAcc_times, 1);
        shadedErrorBar(circle_list, beta_Avg, beta_Sem, {'Color', colorComp(iExp, :), 'MarkerFaceColor', colorComp(iExp, :), 'LineStyle', '-', 'LineWidth', 3}, 0.5); hold on;
        ylim([0, 1]);
        ylimit = ylim;
        %% FDR correction
        % stat_mat = zeros(length(circle_list), 2); % 2: p value and t score
        [~, ~, ~, adj_p] = fdr_bh(stat_mat(:, 1), 0.05, 'pdep'); % method: 'dep', 'pdep'
        for iTime = 1 : length(circle_list)
            p_i = adj_p(iTime);
            if p_i < 0.01
                plot(circle_list(iTime), xLoc * ylimit(end), 'Marker', '*', 'MarkerSize', 4, 'Color', colorComp(iExp, :)); hold on;
            elseif p_i < 0.05 && p_i >= 0.01
                plot(circle_list(iTime), xLoc * ylimit(end), 'Marker', '.', 'MarkerSize', 10, 'Color', colorComp(iExp, :), 'MarkerFaceColor', colorComp(iExp, :), 'LineStyle', 'none', 'LineWidth', 1.5); hold on;
            end
        end
        xLoc = xLoc - 0.05;
    end
    %xlim([0, 0.8]);
    plot(xlim, [(1/2+1/3+1/4)/3, (1/2+1/3+1/4)/3], 'k--', 'LineWidth', 1); hold on;
    plot([0.8, 0.8], ylim, 'k--', 'LineWidth', 1); hold on;
    set(gca, 'LineWidth', 2);
    set(gca, 'FontSize', 18, 'FontWeight', 'bold', 'FontName', 'Arial');
    box off;
    
end



%% accuracy comparison between Exp1.1 (implicit) and Exp1.4 (implicit+cluster)
%%% only the Hamiltonian trials
figure('Position', [100 100 400 300]), clf;
expIds = [1, 2, 3];
for i = 1 : length(expIds)
    angAcc_exp_i = angAcc_Rand_exp(:, :, expIds(i)); % angAcc_Hami_exp(:, :, expIds(i)); % 
    if expIds(i) == 4 
        angAcc_exp_i(12 : end, :) = [];
    elseif expIds(i) == 5
        angAcc_exp_i(14 : end, :) = [];
    end
    [acc_avg, acc_sem] = Mean_and_Se(angAcc_exp_i, 1);
    shadedErrorBar(circle_list, acc_avg, acc_sem, {'Color', colorSets(i, :), 'MarkerFaceColor', colorSets(i, :), 'LineStyle', '-', 'LineWidth', 2}, 0.5); hold on;
    ylim([0, 1]);
end
plot(xlim, [(1/2+1/3+1/4)/3, (1/2+1/3+1/4)/3], 'k--', 'LineWidth', 1); hold on;
plot([0.8, 0.8], ylim, 'k:', 'LineWidth', 1); hold on;
axis xy;
set(gca, 'LineWidth', 2);
set(gca, 'FontSize', 15, 'FontWeight', 'bold', 'FontName', 'Arial');
box off;

%% -------------Part 2: separately for 1, 2 & 3 distractors-------------
%% statistical testing in each time point
stat_mat_dtr = zeros(length(circle_list), 2, 3, length(expList)); % 3 means distractors
for iExp = 1 : length(expList)
    for iTime = 1 : length(circle_list)
        for iDtr = 1 : 3
            angAcc_exp_i = angAcc_dtr_exp(:, iTime, iDtr, iExp); % angAcc_dtr_exp = zeros(subLen, length(circle_list), 3, length(expList));
            [h, p, ci, stats] = ttest(angAcc_exp_i, 1 / (iDtr + 1), 'Tail', 'right');
            stat_mat_dtr(iTime, 1, iDtr, iExp) = p;
            stat_mat_dtr(iTime, 2, iDtr, iExp) = stats.tstat;
        end
    end
end

%% plotting the accuracy time curves
figure('Position', [100 100 1000 300]), clf;
ha = tight_subplot(1, length(expList), [.05 .05], [.1 .02], [.05 .05]);
ha_i = 1;
for iExp = 1 : length(expList)
    axes(ha(ha_i));
    xLoc = 0.98;
    for iDtr = 1 : 3
        angAcc_exp_i = angAcc_dtr_exp(:, :, iDtr, iExp);
        [acc_avg, acc_sem] = Mean_and_Se(angAcc_exp_i, 1);
        
        shadedErrorBar(circle_list, acc_avg, acc_sem, {'Color', colorSets(iDtr, :), 'MarkerFaceColor', colorSets(iDtr, :), 'LineStyle', '-', 'LineWidth', 2}, 0.5); hold on;
        
        %% labeling the significance
        ylim([0, 1]);
        ylimit = ylim;
        [~, ~, ~, adj_p] = fdr_bh(squeeze(stat_mat_dtr(:, 1, iDtr, iExp)), 0.05, 'pdep'); % method: 'dep', 'pdep'
        for iCir = 1 : length(circle_list)
            pval_j = adj_p(iCir);
            if pval_j < 0.01
                plot(circle_list(iCir), xLoc * ylimit(end), 'Marker', '.', 'MarkerSize', 8, 'Color', colorSets(iDtr, :), 'MarkerFaceColor', colorSets(iDtr, :), 'LineStyle', 'none', 'LineWidth', 1.5); hold on;
            elseif pval_j < 0.05 && pval_j >= 0.01
                plot(circle_list(iCir), xLoc * ylimit(end), 'Marker', '>', 'MarkerSize', 4, 'Color', colorSets(iDtr, :), 'MarkerFaceColor', colorSets(iDtr, :), 'LineStyle', 'none', 'LineWidth', 1.5); hold on;
            end
            
        end
        xLoc = xLoc - 0.03;
        plot(xlim, [1/(iDtr + 1), 1/(iDtr + 1)], 'Color', colorSets(iDtr, :), 'LineStyle', '--', 'LineWidth', 1); hold on;
    end
    if iExp == 1 || iExp == 3
        plot([0.8, 0.8], ylim, 'k:', 'LineWidth', 1); hold on;
    end
    
    axis xy;
    set(gca, 'LineWidth', 2);
    set(gca, 'FontSize', 15, 'FontWeight', 'bold', 'FontName', 'Arial');
    box off;
    ha_i = ha_i + 1;
end

%% %% -------------Part 3: random-walk trials (1-700 trials) and hamiltonian-walk trials (701-1500 trials)-------------
% added by rxj @ 03/14/2022
% angAcc_Rand_exp = zeros(subLen, length(circle_list), length(expList)); %% accuracy in each time point
% angAcc_Hami_exp = zeros(subLen, length(circle_list), length(expList)); %% accuracy in each time point
figure('Position', [100 100 1000 300]), clf;
ha = tight_subplot(1, length(expList), [.05 .05], [.1 .02], [.05 .05]);
LineList  = {'-', ':', '-.'};
colorList = {'r', 'b', 'k'};
ha_i = 1;
for iExp = 1 : length(expList)
    axes(ha(ha_i));
    %% plotting the curve for random and hamiltonian walks
    for j = 1 : 3
        if j == 1     %% random
            angAcc_exp_i = angAcc_Rand_exp(:, :, iExp);
        elseif j == 2 %% hamiltonian
            angAcc_exp_i = angAcc_Hami_exp(:, :, iExp);
        elseif j == 3
            angAcc_exp_i = angAcc_exp(:, :, iExp);
        end
        [acc_avg, acc_sem] = Mean_and_Se(angAcc_exp_i, 1);
        shadedErrorBar(circle_list, acc_avg, acc_sem, {'Color', colorList{j}, 'MarkerFaceColor', colorList{j}, 'LineStyle', LineList{j}, 'LineWidth', 2}, 0.5); hold on;
    end
    ylim([0, 1]);
    ylimit = ylim;
    xLoc = 0.98;
    %%
%     [~, ~, ~, adj_p] = fdr_bh(squeeze(stat_mat(:, 1, iExp)), 0.05, 'pdep'); % method: 'dep', 'pdep'
%     for iCir = 1 : length(circle_list)
%         pval_j = adj_p(iCir);
%         if pval_j < 0.01
%             plot(circle_list(iCir), xLoc * ylimit(end), 'Marker', '.', 'MarkerSize', 8, 'Color', colorSets(iExp, :), 'MarkerFaceColor', colorSets(iExp, :), 'LineStyle', 'none', 'LineWidth', 1.5); hold on;
%         elseif pval_j < 0.05 && pval_j >= 0.01
%             plot(circle_list(iCir), xLoc * ylimit(end), 'Marker', '>', 'MarkerSize', 4, 'Color', colorSets(iExp, :), 'MarkerFaceColor', colorSets(iExp, :), 'LineStyle', 'none', 'LineWidth', 1.5); hold on;
%         end
%     end
    if iExp == 1 || iExp == 3
        plot([0.8, 0.8], ylim, 'k:', 'LineWidth', 1); hold on;
    end
    plot(xlim, [(1/2+1/3+1/4)/3, (1/2+1/3+1/4)/3], 'k--', 'LineWidth', 1); hold on;
    axis xy;
    set(gca, 'LineWidth', 2);
    set(gca, 'FontSize', 15, 'FontWeight', 'bold', 'FontName', 'Arial');
    box off;
    
    ha_i = ha_i + 1;
end

%% accuracy for the within- and between-transitions in Random and Hamiltonian Walk trials separately
% write by rxj @ 09/01/2022
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
colorExp = colorSets([1, 3, 5], :);
chance_i = mean([1/2, 1/3, 1/4]);
for iExp = 1 : 3
    colorTmp = colorExp(iExp, :);
    [accAvg_rOh, accSem_rOh] = Mean_and_Se(acc_trans_walk, 1);
    accAvg_rOh = squeeze(accAvg_rOh);
    accSem_rOh = squeeze(accSem_rOh); % acc_trans_walk = zeros(subLen, 2, 2, length(expList));
    figure('Position', [100 100 280 180]), clf;
    barPos = [1, 1.5; 1.7, 2.2];
    for iB = 1 : 2
        barPos_i = barPos(iB, :);
        if iB == 1     %% random walk
            LineBar = '-';
        elseif iB == 2 %% hamiltonian walk
            LineBar = ':';
        end
        %%% line plot
        plot(barPos_i, acc_trans_walk(:, :, iB, iExp), '.-', 'LineWidth', 2, 'MarkerSize', 10, 'Color', [0.6, 0.6, 0.6]); hold on; % individuals
        plot(barPos_i, accAvg_rOh(:, iB, iExp), 'Marker', '.', 'MarkerSize', 10, 'Color', colorTmp, 'LineStyle', 'none'); hold on;
        errorbar(barPos_i(1), accAvg_rOh(1, iB, iExp), accSem_rOh(1, iB, iExp), 'Color', colorTmp, 'Marker', '.', 'LineStyle', '-', 'LineWidth', 4, 'MarkerSize', 10);
        errorbar(barPos_i(2), accAvg_rOh(2, iB, iExp), accSem_rOh(2, iB, iExp), 'Color', colorTmp, 'Marker', '.', 'LineStyle', '-', 'LineWidth', 4, 'MarkerSize', 10);
    end
    set(gca, 'FontSize', 16, 'FontWeight', 'Bold', 'LineWidth', 2);
    set(gca, 'XTick', '', 'XTickLabel', '');
    xlim([0.9, barPos(end, end)+0.1]);
    ylim([0.2, 0.8]);
    plot(xlim, [chance_i, chance_i], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 1); hold on;
    box off;
end

%% paired-samplle t-test for within- and between-transitions
walkNo = 2;
pairT_p = nan(length(expList), walkNo); % walkNo: random and hamiltonian walk
pairT_t = nan(length(expList), walkNo);
for iExp = 1 : length(expList)
    for iB = 1 : 2
        acc_trans_walk_i = acc_trans_walk(:, :, iB, iExp);
        [h, p, ci, stats] = ttest(acc_trans_walk_i(:, 1), acc_trans_walk_i(:, 2));
        pairT_p(iExp, iB) = p;
        pairT_t(iExp, iB) = stats.tstat;
    end
end

%% one-sample ttest: comparison with chance level

%% trial proportions for within-cluster and between-cluster transitions under Random- and Hamiltonian-Walk
% added by rxj @ 05/10/2023
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
colorExp = colorSets([1, 3, 5], :);
figure('Position', [100 100 250 250]), clf;
hold on;
for iExp = 1 : 3
    if iExp == 1
        xLine = [1, 2, 3, 4];
    elseif iExp == 2
        xLine = [1.2, 2.2, 3.2, 4.2];
    elseif iExp == 3
        xLine = [1.4, 2.4, 3.4, 4.4];
    end
    data_trlPro = trlPro_trans_walk(:, :, :, iExp);
    [proAvg, proSem] = Mean_and_Se(data_trlPro, 1);
    proAvg = squeeze(proAvg); % 1st dimension: within- & between-cluster transitions; 2nd dimension: Random & Hamiltonian Walk
    proSem = squeeze(proSem);
    errorbar(xLine(1 : 2), proAvg(:, 1), proSem(:, 1), 'Marker', '.', 'MarkerSize', 20, 'Color', colorExp(iExp, :), 'LineStyle', '-', 'LineWidth', 2); hold on;
    errorbar(xLine(3 : 4), proAvg(:, 2), proSem(:, 2), 'Marker', '.', 'MarkerSize', 20, 'Color', colorExp(iExp, :), 'LineStyle', '-', 'LineWidth', 2); hold on;
    
end
set(gca, 'FontSize', 18, 'FontWeight', 'Bold', 'LineWidth', 2);
set(gca, 'XTick', 1 : 1 : 4, 'XTickLabel', '');
xlim([0, 4.4 + 1]);
box off;





