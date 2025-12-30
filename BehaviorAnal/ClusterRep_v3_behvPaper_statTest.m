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
cd ../
addpath('tight_subplot/');
cd 'HierarchicalCluster/'

addpath('clusterRep-v2-dataAnal/');
addpath('clusterRep-v3-dataAnal/');
addpath('clusterRep-version2/');

%% Subject
ExpWord_List = {'cue', 'explicit', 'cueRand'};
ExpIdx       = 1;
ExpWord      = ExpWord_List{ExpIdx};
if isequal(ExpWord, 'cue')
    subj_list = {'wsn_1_f_18', 'dy_2_f_22', 'haq_3_f_24', 'hry_4_f_20', 'zjx_5_m_20', 'yyq_6_f_18', 'zkw_7_m_18', 'zy_8_f_20', 'hys_9_m_20', 'cjj_10_m_18', ...
                 'dwq_11_f_22', 'ljl_12_m_20', 'jyx_13_m_19', 'zk_14_f_21', 'lsy_15_m_19', 'cjl_16_m_19', 'yjy_17_f_23', 'lym_18_f_19', 'pr_19_f_23', 'ws_20_f_21', ...
                 'wn_21_f_21', 'hjy_22_f_18', 'qyk_23_f_22', 'yd_24_f_20'};
    subjLab   = {'wsn1', 'dy2', 'haq3', 'hry4', 'zjx5', 'yyq6', 'zkw7', 'zy8', 'hys9', 'cjj10', ...
                 'dwq11', 'ljl12', 'jyx13', 'zk14', 'lsy15', 'cjl16', 'yjy17', 'lym18', 'pr19', 'ws20', ...
                 'wn21', 'hjy22', 'qyk23', 'yd24'};
elseif isequal(ExpWord, 'explicit')
    subj_list = {'hsp_1_m_21', 'hr_2_f_22', 'pxy_3_m_19', 'wyx_4_m_22', 'lml_5_f_19', 'lf_6_m_21', 'md_7_f_21', 'srz_8_m_20', 'fsq_9_f_18', 'caq_10_f_20', ...
                 'xjm_11_f_26', 'cjs_12_m_19', 'xxx_13_f_19', 'wsq_14_f_19', 'zzm_15_f_18', 'lsy_16_f_20', 'man_17_f_19', 'zxy_18_f_18', 'cyh_19_f_21', 'szn_20_f_19', ...
                 'sb_21_m_22'}; 
    subjLab   = {'hsp1', 'hr2', 'pxy3', 'wyx4', 'lml5', 'lf6', 'md7', 'srz8', 'fsq9', 'caq10', ...
                 'xjm11', 'cjs12', 'xxx13', 'wsq14', 'zzm15', 'lsy16', 'man17', 'zxy18', 'cyh19', 'szn20', ...
                 'sb21'};
elseif isequal(ExpWord, 'cueRand')
    subj_list = {'zyh_1_f_21', 'why_2_m_18', 'cr_3_f_19', 'zyx_4_f_25', 'wym_5_m_23', 'wd_6_f_18', 'lyh_7_f_21', 'zr_8_f_24', 'zyh_9_m_19', 'zzy_10_f_19', ...
                 'smq_11_f_25', 'sz_12_f_20', 'lzy_13_m_22', 'yxy_14_f_18', 'zxl_15_f_18', 'wqh_16_m_23', 'zxj_17_m_21'}; % , 'xr_18_f_23'
    subjLab   = {'zyh1', 'why2', 'cr3', 'zyx4', 'wym5', 'wd6', 'lyh7', 'zr8', 'lyh9', 'zzy10', ...
                 'smq11', 'sz12', 'lzy13', 'yxy14', 'zxl15', 'wqh16', 'zxj17'}; % , 'xr18'
end
subjNum   = 1 : length(subj_list);
subj_listBv = subj_list;
subLen      = length(subj_list);

%% Experiment parameters
folder     = pwd;
if isequal(ExpWord, 'cue')
    bhvDataDir = [folder, '/FormalExp-LynnNetwork-Results/'];
elseif isequal(ExpWord, 'explicit')
    bhvDataDir = [folder, '/FormalExp-ExplicitLearning-LynnNetwork-Results/'];
elseif isequal(ExpWord, 'cueRand')
    bhvDataDir = [folder, '/FormalExp-LynnNetwork-ImplicitRandom-Results/'];
end
rndTrial = 700;
HamTrial = 800; 
nTrials  = rndTrial + HamTrial;
rndInter = 85;
HamInter = 15;

tInBlc = [230, 2*230, 2*230+240, rndTrial + (200 : 200 : HamTrial)];
nBlock = length(tInBlc); %% seperate the total trials into 10 blocks
trialsInBlc = zeros(nBlock, 2);
trialsInBlc(:, 1) = [1, tInBlc(1 : end - 1) + 1]';
trialsInBlc(:, 2) = tInBlc;

tgtNum  = 1;
dtrNums = [1, 2, 3];
nodeNum = 15;
clsDef  = [1, 2, 3, 4, 5; ...  %% the 1st cluster
           6, 7, 8, 9, 10; ... %% the 2nd cluster
           11, 12, 13, 14, 15];%% the 3rd cluster
transMat = LynnNetwork();

% loss function
marioWin = 100; % For each trial, the maximal bonus is 100 cents(?).
marioLos = -marioWin;
lf_slope = -marioWin/1.5;      %% slope of linear function, y = ax + b
lf_inter = 100 + marioWin/1.5; %% intercept of linear function

% Durations
fixationDur = 0.6; %% seconds
respDur     = 2;   %3s;
grayoffT    = 0.8; %1s;   %% For the 1st second, the gray cue won't appear.
grayOnT     = 0.8; %1s;   %% For 1-2s, the gray square will be overlapped on the target image; for 2-3s, the gray square superimposed on the target image with the largest gray scale
bonusDec    = 0.5; %% bonus decreased from 100 to 0 within 0.5s
bonusTime   = 1;   %% Only for practice trials

%% Columns index for the clusterResult.mat 
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

%% Analyzing the data
lme_subj = cell(1, 2);
for SubIdx = 1 : subLen
    subjBv      = subj_listBv{SubIdx};
    subTit      = subjLab{SubIdx};
    subjDir     = [bhvDataDir, subjBv, '/'];
    load([subjDir, subjBv, 'clusterResult_Blocks.mat'], 'clusterResult');
    trials_Col   = (1 : 1 : nTrials)';
    respRT_Col   = clusterResult(:, respRT);
    errorIdx_Col = clusterResult(:, errorIdx);
    bonus_Col    = clusterResult(:, bonusCn);
    rndHam_Col   = clusterResult(:, trialTp); 
    fast_Col     = clusterResult(:, fastIdx); 
    rndHam_Col(rndHam_Col == 1) = 0; %% random trial
    rndHam_Col(rndHam_Col == 2) = 1; %% hamiltonian trial
    corrIdx      = find(errorIdx_Col ~= 1 & fast_Col ~= 1);
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
    %%% table for the lme test: RT, 
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
        %%% if only the last recency
%         if isempty(rec_iT)
%             recency_nodes(iT) = nodeNum;
%         else 
%             recency_nodes(iT) = iT - rec_iT(end);
%         end
    end
    dtNum = arrayfun(@(x) length(find(dt_nodes(x, :)~=0)), (1 : nTrials)');
    is_BdNode = zeros(nTrials, 1); %% whether the "from" nodes are boundary nodes
    is_BdNode(find(from_nodes == 1 | from_nodes == 5 | from_nodes == 6 | from_nodes == 10 | from_nodes == 11 | from_nodes == 15)) = 1;
    to_BdNode = zeros(nTrials, 1); %% whether the "to" nodes are boundary nodes
    to_BdNode(find(to_nodes == 1 | to_nodes == 5 | to_nodes == 6 | to_nodes == 10 | to_nodes == 11 | to_nodes == 15)) = 1;
    %%% 4 kinds of codings
    lme_h = [];
    lme_h = [repmat(SubIdx, nTrials, 1), trials_Col, respRT_Col, dtNum, is_BdNode, to_BdNode, recency_nodes, rndHam_Col, lureIn];
    lme_h = lme_h(corrIdx, :);
    lme_h(:, 2) = log(lme_h(:, 2)); %% log(trial)
    %%% 3 kinds of codings
    lme_i = [];
    lme_i = [repmat(SubIdx, nTrials, 1), trials_Col, respRT_Col, dtNum, is_BdNode, recency_nodes, transStyle, rndHam_Col, lureIn];
    lme_i = lme_i(corrIdx, :);
    lme_i(:, 2) = log(lme_i(:, 2)); %% log(trial)
    
    lme_subj{1} = [lme_subj{1}; lme_h];
    lme_subj{2} = [lme_subj{2}; lme_i];
end

%% Linear mixed effect model
%% (1) using 'from' and 'to' to code 4 transitions
colName = {'subID', 'logtrialNo', 'rt', 'dtNum', 'fromBdNode', 'toBdNode', 'recencyL', 'recencySL', 'recencyTL', 'rndORham', 'lureWith'};
RT_tbl = array2table(lme_subj{1}, 'VariableNames', colName);
RT_tbl.fromBdNode = nominal(RT_tbl.fromBdNode);
RT_tbl.toBdNode   = nominal(RT_tbl.toBdNode);
RT_tbl.rndORham   = nominal(RT_tbl.rndORham);
RT_tbl.lureWith   = nominal(RT_tbl.lureWith);
% or RT_tbl = array2table(lme_subj); RT_tbl.Properties.VariableNames(1 : 8) = colName;
lmeOutput1 = fitlme(RT_tbl, 'rt~logtrialNo + dtNum + recencyL + recencySL + rndORham + fromBdNode*toBdNode + (1 + logtrialNo + dtNum + recencyL + recencySL + rndORham + fromBdNode*toBdNode|subID)') %  + logtrial*stage + recency + clusteridx
anova(lmeOutput1, 'DFMethod', 'satterthwaite') % 'DFMethod': 'residual' (default), 'satterthwaite', 'none'A

%% (2) using 'isBdNode' and 'transType' to code 3 transitions
colName = {'subID', 'logtrialNo', 'rt', 'dtNum', 'isBdNode', 'recencyL', 'recencySL', 'transType', 'rndORham'};
RT_tbl = array2table(lme_subj{2} , 'VariableNames', colName);
lmeOutput2 = fitlme(RT_tbl, 'rt~logtrialNo + dtNum + isBdNode*transType + recencyL + recencySL + rndORham + (1 + logtrialNo + isBdNode*transType + recencyL + recencySL|subID)'); %  + logtrial*stage + recency + clusteridx

%% calcualte the correlation matrix
[R1, P1] = corrcoef([lme_subj{2}(:, [5, 6]), lme_subj{2}(:, 5) .* lme_subj{2}(:, 6)])
[R2, P2] = corrcoef([lme_subj{2}(:, [5, 8]), lme_subj{2}(:, 5) .* lme_subj{2}(:, 8)])


%% if testing lure distractor effects
lme_lure = lme_subj{1};
lme_lure = lme_lure(find(lme_lure(:, 5) == 1), :);
lme_lure(:, 5) = [];

colName = {'subID', 'logtrialNo', 'rt', 'dtNum', 'toBdNode', 'recencyL', 'recencySL', 'recencyTL', 'rndORham', 'lureWith'};
RT_tbl = array2table(lme_lure, 'VariableNames', colName);
RT_tbl.toBdNode   = nominal(RT_tbl.toBdNode);
RT_tbl.rndORham   = nominal(RT_tbl.rndORham);
RT_tbl.lureWith   = nominal(RT_tbl.lureWith);

lmeOutput1 = fitlme(RT_tbl, 'rt~logtrialNo + dtNum + recencyL + recencySL + rndORham + lureWith*dtNum + lureWith*toBdNode + (1 + logtrialNo + recencyL + recencySL + rndORham + lureWith*dtNum + lureWith*toBdNode|subID)') %  + logtrial*stage + recency + clusteridx
anova(lmeOutput1, 'DFMethod', 'satterthwaite') % 'DFMethod': 'residual' (default), 'satterthwaite', 'none'A



