% ClusterRep_v3_mouseTraj_regressModel
% write by rxj @ 06/15/2022
% using variables estimated from best models as regressors to predict
% subjects' behavior

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

%% Extracting relevant parameters from the model
% added by rxj @ 08/12/2021
fitWord   = 'angleChoice'; % {'clickChoice', 'angleChoice'}
iteWord   = 'naive';       % {'naive', 'trained'}
trialWord = 'all';         % {'all', 'Hamiltonian'}
states = 15;
nFit   = 100;
expList   = {'ImplicitExp', 'ExplicitExp', 'ImplicitRandExp'};
MList     = {'SRheurProcessTraceTruncAttMod'}; 
priorList = {'nullPrior'}; 
folder    = '/Users/ren/Projects-NeuroCode/MyExperiment/HierarchicalCluster';
refit     = 0;
calFlg    = 0; % 0: one-time; 1: accumulated
nanRemove = 1; % 1: remove 'nan' from choiceId and other variales before importing them into fitting procedure
nSim = 500;

%%
circle_list = 0 : 1/60 : 1.5;
colName_bc  = {'Choice', 'trialCnt', 'probTrans', 'probAss', 'EVtrans', 'EVass'}; % binary choice
paramsName_Prob = {'trialCnt', 'probTrans', 'probAss', 'probTrans:probAss'}; 
paramsName_EVs  = {'trialCnt', 'EVtrans', 'EVass', 'EVtrans:EVass'}; 
nParams_Prob    = length(paramsName_Prob); % exclude intercept parameter
nParams_EVs     = length(paramsName_EVs);
subLen = 24;

betas_times_nExp            = cell(subLen, length(expList));
betas_times_firstHalf_nExp  = cell(subLen, length(expList));
betas_times_secondHalf_nExp = cell(subLen, length(expList));
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
        
    elseif isequal(ExpWord, 'ExplicitExp')
        bhvDataDir = [folder, '/FormalExp-ExplicitLearning-LynnNetwork-Results/'];
        subj_list  = {'hsp_1_m_21', 'hr_2_f_22', 'pxy_3_m_19', 'wyx_4_m_22', 'lml_5_f_19', 'lf_6_m_21', 'md_7_f_21', 'srz_8_m_20', 'fsq_9_f_18', 'caq_10_f_20', ...
                      'xjm_11_f_26', 'cjs_12_m_19', 'xxx_13_f_19', 'wsq_14_f_19', 'zzm_15_f_18', 'lsy_16_f_20', 'man_17_f_19', 'zxy_18_f_18', 'cyh_19_f_21', 'szn_20_f_19', ...
                      'sb_21_m_22', 'drq_22_f_18', 'cy_23_f_18', 'zhc_24_m_21'};
        subjLab    = {'hsp1', 'hr2', 'pxy3', 'wyx4', 'lml5', 'lf6', 'md7', 'srz8', 'fsq9', 'caq10', ...
                      'xjm11', 'cjs12', 'xxx13', 'wsq14', 'zzm15', 'lsy16', 'man17', 'zxy18', 'cyh19', 'szn20', ...
                      'sb21', 'drq22', 'cy23', 'zhc24'};
        
    elseif isequal(ExpWord, 'ImplicitRandExp')
        bhvDataDir = [folder, '/FormalExp-LynnNetwork-ImplicitRandom-Results/'];
        subj_list  = {'zyh_1_f_21', 'why_2_m_18', 'cr_3_f_19', 'zyx_4_f_25', 'wym_5_m_23', 'wd_6_f_18', 'lyh_7_f_21', 'zr_8_f_24', 'zyh_9_m_19', 'zzy_10_f_19', ...
                      'smq_11_f_25', 'sz_12_f_20', 'lzy_13_m_22', 'yxy_14_f_18', 'zxl_15_f_18', 'wqh_16_m_23', 'zxj_17_m_21', 'skx_18_f_20', 'zlh_19_m_24', 'gwt_20_f_23', ...
                      'lwn_21_f_18', 'lrp_22_m_21', 'sjj_23_f_18', 'xy_24_f_19'}; % , 'xr_18_f_23'
        subjLab    = {'zyh1', 'why2', 'cr3', 'zyx4', 'wym5', 'wd6', 'lyh7', 'zr8', 'lyh9', 'zzy10', ...
                      'smq11', 'sz12', 'lzy13', 'yxy14', 'zxl15', 'wqh16', 'zxj17', 'skx18', 'zlh19', 'gwt20', ...
                      'lwn21', 'lrp22', 'sjj23', 'xy24'}; % , 'xr18'

    end
    tInBlc = [230, 2*230, 2*230+240, rndTrial + (200 : 200 : HamTrial)];
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
    
    %% Loop over participants

    for iSub = 1 : subLen
        disp([ExpWord, '-subj', num2str(iSub)]);
        %%
        subID  = subjLab{iSub};
        subjBv = subj_listBv{iSub};
        subjDir = [bhvDataDir, subjBv, '/'];
        
        if isequal(expMode_i, 'mouse')
            load([subjDir, subjBv, 'clusterResult_Blocks.mat'], 'clusterResult');
 
            trials_Col   = (1 : 1 : nTrials)';
            respRT_Col   = clusterResult(:, respRT);
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

            load([subjDir, subjBv,  'respYes_trials_blc.mat'], 'respYes_trials_blc');
            load([subjDir, subjBv,  'mouseTraj_trials_blc.mat'], 'mouseTraj_trials_blc');

            %%% extract the angles
            choiceId_ang  = nan(length(stim), 1);
            choiceId_time = nan(length(stim), length(circle_list));
            dtrNoCnt_iSub = nan(length(stim), 1);
            objAng_col    = zeros(length(stim), 4);
            objDtrNo_col  = zeros(length(stim), 1);
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
                    objAng_col(iCount, :) = objAngs_stay;
                    objDtrNo_col(iCount)  = dtrNoCnt_j(1);

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
                                    choiceId_time(iCount, iTp) = 1;
                                    
                                elseif minId ~= 1
                                    choiceId_time(iCount, iTp) = 0;
                                    
                                end
                            end
                        end
                    end
                    iCount = iCount + 1;
                end
            end
            
            %% classify 'Yes' or 'No' according to angles in mouse trajectory or subjects' selections
            choiceId = choiceId_ang;
            if nanRemove == 1
                trials_Col(isnan(choiceId)) = [];
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
                choiceId_time(isnan(choiceId), :)     = [];
                choiceId(isnan(choiceId))             = [];
            end

            %% fitting procedure
            nData_i = length(stim);
            if isequal(iteWord, 'naive')       %% learn the transition matrix from zero
                Aest_input = eye(states, states);
                nij_input  = zeros(states, states);
                M_input    = eye(states, states);
            end

            % ------ Only one model ------
            Midx      = MList{1};
            priorWord = priorList{1}; %'InformPrior', 'nullPrior';

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

            [~, paramsEst] = ClusterRep_FitModels_choice(ExpWord, subID, stim, resp, states, choiceId, delta_t_trials, assDist_mat, fitWord, ...
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
                [~, M_trans, M_ass, p_choice_all, choiceID_all, lambda_trials, p_choice_trs, p_choice_ass, lr_trials, EV_trs, EV_ass] = ...
                    SR_heurProcessTrace_choice_gdChoice_truncAttMod(paramsEst, stim, resp, states, choiceId, delta_t_trials, assDist_mat, expMode_i, posXY, respYes_trials_blc, ...
                                                        mouseTraj_trials_blc, objAng_col, objDtrNo_col, miniDmat, fitWord, transMat_input, nij_input, priorWord);
            end

            %% ---------- Model-based regression analysis ----------
            %% All trials
            % 1) whether the prediction-based Prob and memory-based Prob could
            % be dissociated across time points
            % 2) EV as the regressor, though the scale predicted from the
            % two system might be different
            %regressor_trial = [log(trials_Col), p_choice_trs(:, 1), p_choice_ass(:, 1), zscore(EV_trs(:, 1)), zscore(EV_ass(:, 1))]; 
            regressor_trial_raw = [log(trials_Col), p_choice_trs(:, 1), p_choice_ass(:, 1), EV_trs(:, 1), EV_ass(:, 1)]; 
            regressor_trial     = regressor_trial_raw;
            regressor_trial(:, 2 : end) = zscore(regressor_trial(:, 2 : end));

            betas_iTime = nan(nParams_Prob, length(circle_list), 2); % 2: probability and EVs as regressors
            for iTime = 1 : length(circle_list)
                choiceId_i  = choiceId_time(:, iTime);
                % ------ Constructing the GLM array ------
                logit_trial = [choiceId_i, regressor_trial];
                logit_trial(isnan(logit_trial(:, 1)), :) = [];
                
                if size(logit_trial, 1) > (3 * nParams_Prob)
                    %%% convert matrix into table
                    lm_tbl = array2table(logit_trial, 'VariableNames', colName_bc);

                    %%% ---------- 1. Regression with predicted probability as the regressors ----------
                    modelspec = 'Choice ~ trialCnt + probTrans + probAss + probTrans*probAss';
                    opts = statset('glmfit');
                    opts.MaxIter = nParams_Prob*1000;
                    fit = fitglm(lm_tbl, modelspec, 'Distribution', 'binomial', 'link', 'logit', 'options', opts); % 
                    fitNames = fit.CoefficientNames;
                    betas_Prob = nan(nParams_Prob, 1);
                    for iPs = 1 : nParams_Prob
                        params_iPs = paramsName_Prob{iPs};
                        if any(strcmp(fitNames, params_iPs))
                            betas_Prob(iPs) = fit.Coefficients{params_iPs, 1}; % 1 is the column of the Estimate
                            % betas = fit.Coefficients.Estimate;
                        end
                    end
                    betas_iTime(:, iTime, 1) = betas_Prob;

                    %%% ---------- 2. Regression with Expected Values as the regressors ----------
                    modelspec = 'Choice ~ trialCnt + EVtrans + EVass + EVtrans*EVass';
                    opts = statset('glmfit');
                    opts.MaxIter = nParams_EVs*1000;
                    fit = fitglm(lm_tbl, modelspec, 'Distribution', 'binomial', 'link', 'logit', 'options', opts); % 
                    fitNames = fit.CoefficientNames;
                    betas_EVs = nan(nParams_EVs, 1);
                    for iPs = 1 : nParams_EVs
                        params_iPs = paramsName_EVs{iPs};
                        if any(strcmp(fitNames, params_iPs))
                            betas_EVs(iPs) = fit.Coefficients{params_iPs, 1}; % 1 is the column of the Estimate
                            % betas = fit.Coefficients.Estimate;
                        end
                    end
                    betas_iTime(:, iTime, 2) = betas_EVs;

                end
            end
            betas_times_nExp{iSub, iExp} = betas_iTime;

            %% ------ First 700 trials (only random walks) ------
            trainTrl_idx = find(trials_Col <= 700);
            regressor_train_trial = regressor_trial_raw(trainTrl_idx, :);
            regressor_train_trial(:, 2 : end) = zscore(regressor_train_trial(:, 2 : end));

            betas_iTime_train = nan(nParams_Prob, length(circle_list), 2); % 2: probability and EVs as regressors
            for iTime = 1 : length(circle_list)
                choiceId_i  = choiceId_time(:, iTime);
                choiceId_i  = choiceId_i(trainTrl_idx);
                % ------ Constructing the GLM array ------
                logit_trial = [choiceId_i, regressor_train_trial];
                logit_trial(isnan(logit_trial(:, 1)), :) = [];
                
                if size(logit_trial, 1) > (3 * nParams_Prob)
                    %%% convert matrix into table
                    lm_tbl = array2table(logit_trial, 'VariableNames', colName_bc);

                    %%% ---------- 1. Regression with predicted probability as the regressors ----------
                    modelspec = 'Choice ~ trialCnt + probTrans + probAss + probTrans*probAss';
                    opts = statset('glmfit');
                    opts.MaxIter = nParams_Prob*1000;
                    fit = fitglm(lm_tbl, modelspec, 'Distribution', 'binomial', 'link', 'logit', 'options', opts); % 
                    fitNames = fit.CoefficientNames;
                    betas_Prob = nan(nParams_Prob, 1);
                    for iPs = 1 : nParams_Prob
                        params_iPs = paramsName_Prob{iPs};
                        if any(strcmp(fitNames, params_iPs))
                            betas_Prob(iPs) = fit.Coefficients{params_iPs, 1}; % 1 is the column of the Estimate
                            % betas = fit.Coefficients.Estimate;
                        end
                    end
                    betas_iTime_train(:, iTime, 1) = betas_Prob;

                    %%% ---------- 2. Regression with Expected Values as the regressors ----------
                    modelspec = 'Choice ~ trialCnt + EVtrans + EVass + EVtrans*EVass';
                    opts = statset('glmfit');
                    opts.MaxIter = nParams_EVs*1000;
                    fit = fitglm(lm_tbl, modelspec, 'Distribution', 'binomial', 'link', 'logit', 'options', opts); % 
                    fitNames = fit.CoefficientNames;
                    betas_EVs = nan(nParams_EVs, 1);
                    for iPs = 1 : nParams_EVs
                        params_iPs = paramsName_EVs{iPs};
                        if any(strcmp(fitNames, params_iPs))
                            betas_EVs(iPs) = fit.Coefficients{params_iPs, 1}; % 1 is the column of the Estimate
                            % betas = fit.Coefficients.Estimate;
                        end
                    end
                    betas_iTime_train(:, iTime, 2) = betas_EVs;

                end
            end
            betas_times_firstHalf_nExp{iSub, iExp} = betas_iTime_train;


            %% ------ Second 800 trials (mixture of random and hamiltonian walks) ------
            testTrl_idx = find(trials_Col > 700);
            regressor_test_trial = regressor_trial_raw(testTrl_idx, :);
            regressor_test_trial(:, 2 : end) = zscore(regressor_test_trial(:, 2 : end));

            betas_iTime_test = nan(nParams_Prob, length(circle_list), 2); % 2: probability and EVs as regressors
            for iTime = 1 : length(circle_list)
                choiceId_i  = choiceId_time(:, iTime);
                choiceId_i  = choiceId_i(testTrl_idx);
                % ------ Constructing the GLM array ------
                logit_trial = [choiceId_i, regressor_test_trial];
                logit_trial(isnan(logit_trial(:, 1)), :) = [];
                
                if size(logit_trial, 1) > (3 * nParams_Prob)
                    %%% convert matrix into table
                    lm_tbl = array2table(logit_trial, 'VariableNames', colName_bc);

                    %%% ---------- 1. Regression with predicted probability as the regressors ----------
                    modelspec = 'Choice ~ trialCnt + probTrans + probAss + probTrans*probAss';
                    opts = statset('glmfit');
                    opts.MaxIter = nParams_Prob*1000;
                    fit = fitglm(lm_tbl, modelspec, 'Distribution', 'binomial', 'link', 'logit', 'options', opts); % 
                    fitNames = fit.CoefficientNames;
                    betas_Prob = nan(nParams_Prob, 1);
                    for iPs = 1 : nParams_Prob
                        params_iPs = paramsName_Prob{iPs};
                        if any(strcmp(fitNames, params_iPs))
                            betas_Prob(iPs) = fit.Coefficients{params_iPs, 1}; % 1 is the column of the Estimate
                            % betas = fit.Coefficients.Estimate;
                        end
                    end
                    betas_iTime_test(:, iTime, 1) = betas_Prob;

                    %%% ---------- 2. Regression with Expected Values as the regressors ----------
                    modelspec = 'Choice ~ trialCnt + EVtrans + EVass + EVtrans*EVass';
                    opts = statset('glmfit');
                    opts.MaxIter = nParams_EVs*1000;
                    fit = fitglm(lm_tbl, modelspec, 'Distribution', 'binomial', 'link', 'logit', 'options', opts); % 
                    fitNames = fit.CoefficientNames;
                    betas_EVs = nan(nParams_EVs, 1);
                    for iPs = 1 : nParams_EVs
                        params_iPs = paramsName_EVs{iPs};
                        if any(strcmp(fitNames, params_iPs))
                            betas_EVs(iPs) = fit.Coefficients{params_iPs, 1}; % 1 is the column of the Estimate
                            % betas = fit.Coefficients.Estimate;
                        end
                    end
                    betas_iTime_test(:, iTime, 2) = betas_EVs;

                end
            end
            betas_times_secondHalf_nExp{iSub, iExp} = betas_iTime_test;
            
            %%
            
        elseif isequal(expMode_i, 'key')
            % ......
 
        end
    end
end

%% Extracting the relevant betas
% betas_times_nExp            = cell(subLen, length(expList));
% betas_times_firstHalf_nExp  = cell(subLen, length(expList));
% betas_times_secondHalf_nExp = cell(subLen, length(expList));

trialFlg = 3; % 1-all trials; 2-first 700 trials; 3-last 800 trials;
if trialFlg == 1
    betas_times_tmp = betas_times_nExp;
elseif trialFlg == 2
    betas_times_tmp = betas_times_firstHalf_nExp;
elseif trialFlg == 3
    betas_times_tmp = betas_times_secondHalf_nExp;
end

dataFlg  = 2; % 1-predicted choice probability as the regressor; 2-EVs
betas_times_plot = nan(3, length(circle_list), subLen, length(expList)); % 3: betas from memory and prediction system and their interactions
% paramsName_Prob = {'trialCnt', 'probTrans', 'probAss', 'probTrans:probAss'}; 
% paramsName_EVs  = {'trialCnt', 'EVtrans', 'EVass', 'EVtrans:EVass'}; 
for iExp = 1 : length(expList)
    for iSub = 1 : subLen
        betas_times_plot(:, :, iSub, iExp) = betas_times_tmp{iSub, iExp}(2 : end, :, dataFlg);
    end
end

%% plotting the beta time curves
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

colorExp = [251, 154, 153; ...
            178, 223, 138; ...
            166, 206, 227] ./ 255;
LineSty = '-';
for iExp = 1 : length(expList)
    disp(['---------- ', expList{iExp}, ' ----------']);
    figure('Position', [100 100 300 140]), clf;

    betas_times_ii = betas_times_plot(:, :, :, iExp);
    [acc_avg, acc_sem] = Mean_and_Se(betas_times_ii, 3); % (3, length(circle_list))

    for ii = 1 : 3
        shadedErrorBar(circle_list, acc_avg(ii, :), acc_sem(ii, :), {'Color', colorExp(ii, :), 'MarkerFaceColor', colorExp(ii, :), 'LineStyle', LineSty, 'LineWidth', 3}, 0.5); hold on;
    end
    ylim([-0.2, 0.8]);
    plot(xlim, [0, 0], 'k--', 'LineWidth', 0.8); hold on;
    plot([0.8, 0.8], ylim, 'k--', 'LineWidth', 1); hold on;
    axis xy;
    set(gca, 'LineWidth', 0.8);
    set(gca, 'FontSize', 10, 'FontWeight', 'bold', 'FontName', 'Arial');
    set(gca, 'XTick', 0 : 0.4 : 1.5, 'XTickLabel', '');
    set(gca, 'YTick', -0.2 : 0.2 : 0.8, 'YTickLabel', '');
    box off;
    ax = gca;
    save_name = ['Exp', num2str(iExp), '-', expList{iExp}, '-BetaCuves.png'];
    exportgraphics(ax, save_name, 'Resolution', 600);
end


