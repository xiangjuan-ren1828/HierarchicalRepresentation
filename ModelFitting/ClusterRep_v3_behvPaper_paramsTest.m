% ClusterRep_v3_behvPaper_paramsTest.m
% modified by XR based on the followings 
% Check the time scale of the prediction-based and memory-based system
% ==================================================
% ClusterRep_model_paramsTest
% write by rxj @ 05/06/2021 based on the followings
% ==================================================
% LynnNC_paramsTest
% write by rxj @ 07/28/2020

clear
clc

%%
addpath('tight_subplot/');
addpath(genpath('HierarchicalCluster/'));
%addpath(genpath('spm12/'));

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
miniDmat = modelRDM_minimalDist();

%% model comparision
% added by rxj @ 08/12/2021
fitWord   = 'angleChoice'; % {'clickChoice', 'angleChoice'}
iteWord   = 'naive';     % {'naive', 'trained'}
trialWord = 'all'; % {'all', 'Hamiltonian'}
states = 15;
nFit   = 100;
expList         = {'ImplicitRandExp'}; %{'ImplicitExp', 'ExplicitExp', 'ImplicitRandExp', 'keypress', 'keypress_EEG'};
MList           = {'objTP', 'MaxEntropy', 'MaxEntropy', 'SR', 'heurProcessProb', 'heurProcessTraceWin', 'heurProcessTrace', 'SRheurProcessProb', 'SRheurProcessTraceWin', 'SRheurProcessTrace'}; 
priorList       = {'nullPrior', 'nullPrior', 'InformPrior', 'nullPrior', 'nullPrior', 'nullPrior', 'nullPrior', 'nullPrior', 'nullPrior', 'nullPrior'};
folder          = '/Users/cdlab/Nextcloud/MATLAB/MyExperiment/HierarchicalCluster'; %pwd;
%folder          = '/Volumes/HierarchicalCluster';
refit           = 0;
%%
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
        
    elseif isequal(ExpWord, 'keypress')
        bhvDataDir = [folder, '/FormalExp-LynnNetwork-KeyPressResults/'];
        subj_list = {'lkx_2_f_18', 'gwz_4_m_21', 'nwq_5_f_23', 'jsp_6_f_18', 'kly_8_m_19', 'lyy_9_f_20', 'nyr_11_f_18', 'zty_13_m_19', ... 'sch_10_m_20', 
                     'dyc_11_f_20', 'ymy_12_f_23', 'dyj_13_f_20', 'zyz_14_m_19', 'zcy_15_m_25', 'jq_16_f_23', 'xly_17_f_19', 'xhq_18_f_19', 'xjf_19_m_20', 'xyq_20_m_20', ...
                     'hqy_21_m_23', 'lyq_22_f_24', 'lmy_23_f_24', 'fy_24_m_25', 'phq_25_f_21', 'yjq_27_f_23', 'wlx_29_m_20', ... 'xj_28_f_19', 
                     'hry_31_f_20', 'ml_32_m_20', 'zt_33_m_18', 'dl_34_m_23', 'ghf_35_m_19', 'sp_36_m_20', 'yjy_37_m_25', 'yxy_38_f_19'};   
        subjLab   = {'lkx2', 'gwz4', 'nwq5', 'jsp6', 'kly8', 'lyy9',  'nyr11', 'zty13', ... 'sch10',
                     'dyc11', 'ymy12', 'dyj13', 'zyz14', 'zcy15', 'jq16', 'xly17', 'xhq18', 'xjf19', 'xyq20', ...
                     'hqy21', 'lyq22', 'lmy23', 'fy24', 'phq25', 'yjq27', 'wlx29', ... 'xj28', 
                     'hry31', 'ml32', 'zt33', 'dl34', 'ghf35', 'sp36', 'yjy37', 'yxy38'}; 
        nBlock = 7;
                 
    elseif isequal(ExpWord, 'keypress_EEG')
        bhvDataDir = [folder, '/EEG_EyeTrack_KeyPress_FormalExp-LynnNetwork-KeyPressResults/'];
        subj_list = {'lkx_2_f_18', 'gwz_4_m_21', 'nwq_5_f_23', 'jsp_6_f_18', 'kly_8_m_19', 'lyy_9_f_20', 'nyr_11_f_18', 'zty_13_m_19', ... 'sch_10_m_20', 
                     'dyc_11_f_20', 'ymy_12_f_23', 'dyj_13_f_20', 'zyz_14_m_19', 'zcy_15_m_25', 'jq_16_f_23', 'xly_17_f_19', 'xhq_18_f_19', 'xjf_19_m_20', 'xyq_20_m_20', ...
                     'hqy_21_m_23', 'lyq_22_f_24', 'lmy_23_f_24', 'fy_24_m_25', 'phq_25_f_21', 'yjq_27_f_23', 'wlx_29_m_20', ... 'xj_28_f_19', 
                     'hry_31_f_20', 'ml_32_m_20', 'zt_33_m_18', 'dl_34_m_23', 'ghf_35_m_19', 'sp_36_m_20', 'yjy_37_m_25', 'yxy_38_f_19'};   
        subjLab   = {'lkx2', 'gwz4', 'nwq5', 'jsp6', 'kly8', 'lyy9',  'nyr11', 'zty13', ... 'sch10',
                     'dyc11', 'ymy12', 'dyj13', 'zyz14', 'zcy15', 'jq16', 'xly17', 'xhq18', 'xjf19', 'xyq20', ...
                     'hqy21', 'lyq22', 'lmy23', 'fy24', 'phq25', 'yjq27', 'wlx29', ... 'xj28', 
                     'hry31', 'ml32', 'zt33', 'dl34', 'ghf35', 'sp36', 'yjy37', 'yxy38'}; 
        nBlock = 6;
        
    elseif isequal(ExpWord, 'keypress_fullExp')
        
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
    params_allModels = cell(subLen, length(MList));
    nData_allModels  = nan(subLen, length(MList));
    for iSub = 1 : subLen
        %%
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
                
                trl_counts(isnan(choiceId)) = [];
                stim(isnan(choiceId))       = [];
                resp(isnan(choiceId), :)    = [];
                objAng_col(isnan(choiceId), :)   = [];
                objDtrNo_col(isnan(choiceId), :) = [];
                %%% In the incorrect version, choiceId(isnan(choiceId)) =
                %%% [] was run before these following 4 variables
                delta_t_trials(isnan(choiceId), :)    = [];
                posXY(isnan(choiceId), :, :)          = [];
                respYes_trials_blc(isnan(choiceId))   = [];
                mouseTraj_trials_blc(isnan(choiceId)) = [];
                choiceId(isnan(choiceId))             = [];

            end
            
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
                [~, paramsEst] = ClusterRep_FitModels_choice(ExpWord, subID, stim, resp, states, choiceId, delta_t_trials, fitWord, ...
                                                          expMode_i, posXY, respYes_trials_blc, mouseTraj_trials_blc, objAng_col, objDtrNo_col, miniDmat, nFit, Midx, folder, refit, transMat_input, nij_input, iteWord, trialWord, priorWord); 
                params_allModels{iSub, iM} = paramsEst;
                nData_allModels(iSub, iM)  = nData_i;
            end
            
            %%
            
        elseif isequal(expMode_i, 'key')
            % ......
          
        end
    end 
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
colorStats = colorSets([1, 3, 5], :);
colorDiv   = colorSets([1, 3, 5, 2, 4, 6], :);

%% parameters
% MList     = {'objTP', 'MaxEntropy', 'MaxEntropy', 'SR', 'heurProcessProb', 'heurProcessTraceWin', 'heurProcessTrace', 'SRheurProcessProb', 'SRheurProcessTraceWin', 'SRheurProcessTrace'}; 
% priorList = {'nullPrior', 'nullPrior', 'InformPrior', 'nullPrior', 'nullPrior', 'nullPrior', 'nullPrior', 'nullPrior', 'nullPrior', 'nullPrior'};
iM = 10; %input('which model, 1-MaxEntropy, 2-SR, 3-SRgd, 5-SRgdClusterTrace, 6-SRgdClusterWin: '); % 4-SRcw, 
Midx = MList{iM};
priorWord = priorList{iM}; %'InformPrior', 'nullPrior';
%%% transition models
if isequal(Midx, 'objTP')
    paramsName = {'\beta_C'};
elseif isequal(Midx, 'MaxEntropy') && isequal(priorWord, 'nullPrior')
    paramsName = {'\beta_m', '\beta_C'};
elseif isequal(Midx, 'MaxEntropy') && isequal(priorWord, 'InformPrior')
    paramsName = {'\beta_m', '\beta_C', 'initV'};
elseif isequal(Midx, 'SR') 
    paramsName = {'\gamma', '\alpha', '\beta_C'};
%%% association models
elseif isequal(Midx, 'heurProcessProb') 
    paramsName = {'winLen', '\beta_C', '\alpha_{mw}'};
elseif isequal(Midx, 'heurProcessTraceWin')   
    paramsName = {'winLen', '\beta_C', '\alpha_{mw}'};
elseif isequal(Midx, 'heurProcessTrace')  
    paramsName = {'\beta_m', '\beta_C', '\alpha_{mw}'};
%%% mixture of transition and association models
elseif isequal(Midx, 'SRheurProcessProb') 
    paramsName = {'\gamma_{tr}', '\alpha_{tr}', 'winLen', '\alpha_{mw}', '\alpha_{gd}', '\beta_{Ctr}', '\beta_{Cas}'};
elseif isequal(Midx, 'SRheurProcessTraceWin') 
    paramsName = {'\gamma_{tr}', '\alpha_{tr}', 'winLen', '\alpha_{mw}', '\alpha_{gd}', '\beta_{Ctr}', '\beta_{Cas}'};
elseif isequal(Midx, 'SRheurProcessTrace')
    paramsName = {'\gamma_{tr}', '\alpha_{tr}', '\beta_m', '\alpha_{mw}', '\alpha_{gd}', '\beta_{Ctr}', '\beta_{Cas}'};
end
params_test = params_allModels(:, iM);
params_test = cell2mat(params_test');

%% plotting: boxplot
for iP = 1 : length(paramsName)
    pTit = paramsName{iP};
    params_iP = (params_test(iP, :))';
    if isequal(pTit, '\beta_{Ctr}') || isequal(pTit, '\beta_{Cas}')
        params_iP = log(params_iP);
    elseif isequal(pTit, '\beta_m') && isequal(Midx, 'SRheurProcessTrace') && isequal(ExpWord, 'ImplicitExp')
        params_iP([7,18]) = [];
        params_iP = log(params_iP);
    elseif isequal(pTit, '\beta_m') && isequal(Midx, 'SRheurProcessTrace') && isequal(ExpWord, 'ImplicitRandExp')    
        params_iP = log(params_iP);
    end
    
    figure('Position', [100 100 150 200]), clf;
    hold on; 
    Bi = 1;
    xRand = linspace(Bi - 0.3, Bi + 0.3, length(params_iP));
    plot(xRand, params_iP, 'Marker', '.', 'MarkerSize', 30, 'Color', 'k', 'LineStyle', 'none'); hold on;
    
    boxplot(params_iP, 'BoxStyle', 'outline', 'Colors', 'k', 'MedianStyle', 'line', 'Notch', 'off', 'symbol', '', 'Whisker', 1.5); hold on; % whisker = 1.5 is the default
    hh = findobj(gca,'Tag','Box');
    for Bi = 1 : length(hh)
        patch(get(hh(Bi),'XData'),get(hh(Bi),'YData'), [0.6, 0.6, 0.6], 'FaceAlpha', 0.6);
    end
    set(findobj(gca, 'type', 'line'), 'LineWidth', 3);
    lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
    set(lines, 'Color', 'k', 'LineWidth', 6);
    xlim([0.5, 1.5]);
    set(gca, 'LineWidth', 2, 'FontSize', 18, 'FontWeight', 'Bold');
    set(gca, 'XTick', 1, 'Box', 'off');
    set(gca, 'XTickLabel','');
    %title(['log(', pTit, ')'], 'FontSize', 20);
    title(pTit, 'FontSize', 20);
end

%% plotting: distribution
for iP = 1 : length(paramsName)
    pTit = paramsName{iP};
    params_iP = (params_test(iP, :))';
    if isequal(pTit, '\beta_{Ctr}') || isequal(pTit, '\beta_{Cas}')
        params_iP = log(params_iP);
    elseif isequal(pTit, '\beta_m') && isequal(Midx, 'SRheurProcessTrace') && isequal(ExpWord, 'ImplicitExp')
        params_iP([7,18]) = [];
        params_iP = log(params_iP);
    elseif isequal(pTit, '\beta_m') && isequal(Midx, 'SRheurProcessTrace') && isequal(ExpWord, 'ImplicitRandExp')    
        params_iP = log(params_iP);
    end
    x = linspace(min(params_iP), max(params_iP), 1000);
    
    figure('Position', [100 100 250 200]), clf;
    hold on; 
    pd_kernel  = fitdist(params_iP, 'Kernel');
    pdf_iP = pdf(pd_kernel, x);
    plot(x, pdf_iP, 'Color', 'k', 'LineWidth', 2); hold on;
    area(x, pdf_iP, 'FaceColor', colorSets(6, :));
    
    xlim([min(params_iP), max(params_iP)]);
    set(gca, 'LineWidth', 2);
    set(gca, 'FontSize', 15, 'FontWeight', 'bold', 'FontName', 'Arial');
    box off;
    title(pTit, 'FontSize', 20);
end

%% Correlation between estimated parameters and subjects' behaviour
% or correlation between estimated parameters and beta coefficients of regression
load([angleDir, 'acc_true_subj.mat'], 'acc_true_subj');
[r, p] = corr(acc_true_subj', (params_test(3, :))')

figure('Position', [100 100 200 200]), clf;
hold on;
plot(acc_true_subj', (params_test(3, :))', 'Marker', '.', 'MarkerSize', 15, 'MarkerFaceColor', [0.6, 0.6, 0.6], 'MarkerEdgeColor', 'k', 'LineStyle', 'none'); hold on;
set(gca, 'LineWidth', 2);
set(gca, 'FontSize', 18, 'FontWeight', 'bold', 'FontName', 'Arial');
box off;





