% ClusterRep_v3_RTAnal
% write by rxj @ 11/03/2020
% =========================================================================
% ClusterRep_v2_RTAnal
% write by rxj @ 09/23/2020
% All of the analysis here are based on the collected data in
% ClusterRep_FormalExp_v4.m

%%
clear
clc

%%
addpath('tight_subplot/');
addpath(genpath('HierarchicalCluster/'));

%% Subject
ExpWord_List = {'cue', 'explicit', 'cueRand', 'cueEyeTrack'};
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
                 'sb_21_m_22', 'drq_22_f_18', 'cy_23_f_18', 'zhc_24_m_21'}; 
    subjLab   = {'hsp1', 'hr2', 'pxy3', 'wyx4', 'lml5', 'lf6', 'md7', 'srz8', 'fsq9', 'caq10', ...
                 'xjm11', 'cjs12', 'xxx13', 'wsq14', 'zzm15', 'lsy16', 'man17', 'zxy18', 'cyh19', 'szn20', ...
                 'sb21', 'drq22', 'cy23', 'zhc24'};  
elseif isequal(ExpWord, 'cueRand')
    subj_list = {'zyh_1_f_21', 'why_2_m_18', 'cr_3_f_19', 'zyx_4_f_25', 'wym_5_m_23', 'wd_6_f_18', 'lyh_7_f_21', 'zr_8_f_24', 'zyh_9_m_19', 'zzy_10_f_19', ...
                 'smq_11_f_25', 'sz_12_f_20', 'lzy_13_m_22', 'yxy_14_f_18', 'zxl_15_f_18', 'wqh_16_m_23', 'zxj_17_m_21', 'skx_18_f_20', 'zlh_19_m_24', 'gwt_20_f_23', ...
                 'lwn_21_f_18', 'lrp_22_m_21', 'sjj_23_f_18', 'xy_24_f_19'}; % , 'xr_18_f_23'
    subjLab   = {'zyh1', 'why2', 'cr3', 'zyx4', 'wym5', 'wd6', 'lyh7', 'zr8', 'lyh9', 'zzy10', ...
                 'smq11', 'sz12', 'lzy13', 'yxy14', 'zxl15', 'wqh16', 'zxj17', 'skx18', 'zlh19', 'gwt20', ...
                 'lwn21', 'lrp22', 'sjj23', 'xy24'}; % , 'xr18'
elseif isequal(ExpWord, 'cueEyeTrack')
    subj_list = {'mcaj_1_f_20', 'yjw_2_f_21', 'lyc_3_m_20', 'cjq_4_f_18', 'lmh_5_f_19', 'lh_6_f_19', 'stj_7_f_20', 'ljl_8_m_21', 'csx_10_m_20', 'ljj_11_f_24', 'lx_12_m_19',  ...
                 'lyc_13_f_17', 'jzc_14_m_24', 'lzz_15_m_22', 'wzx_16_m_22', 'mym_17_f_18', 'ys_18_f_18', 'djh_19_m_18'}; % without eye data: 'mcaj_1_f_20', 'cjq_4_f_18'
    subjLab   = {'mcaj1', 'yjw2', 'lyc3', 'cjq4', 'lmh5', 'lh6', 'stj7', 'ljl8', 'csx10', 'ljj11', 'lx12', ...
                 'lyc13', 'jzc14', 'lzz15', 'wzx16', 'mym17', 'ys18', 'djh19'}; % 'mcaj1', 'cjq4'    
end
subjNum   = 1 : length(subj_list);
subj_listBv = subj_list;
subLen      = length(subj_list);

%% Experiment parameters
%folder     = pwd;
folder = '/Users/renxiangjuan/Nextcloud/MATLAB/MyExperiment/HierarchicalCluster';
if isequal(ExpWord, 'cue')
    bhvDataDir = [folder, '/FormalExp-LynnNetwork-Results/']; 
elseif isequal(ExpWord, 'explicit')
    bhvDataDir = [folder, '/FormalExp-ExplicitLearning-LynnNetwork-Results/'];
elseif isequal(ExpWord, 'cueRand')
    bhvDataDir = [folder, '/FormalExp-LynnNetwork-ImplicitRandom-Results/'];
elseif isequal(ExpWord, 'cueEyeTrack')
    bhvDataDir = [folder, '/EyeTrackFormalExp-LynnNetwork-Results/'];
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
predRatio = zeros(subLen, nBlock);
errRatio  = zeros(subLen, nBlock); %% ratio of error or time out trials
fastRatio = zeros(subLen, nBlock);
corrPreRt = zeros(subLen, 1); %% correction ratio of prediction trials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
baseAcc      = zeros(subLen, 1);
transTpAcc   = zeros(subLen, 2); %% correct response ratio in within- and between-cluster transitions
transTpAcc_lure = zeros(subLen, 2, 2); % the 1st and 2nd '2' denote 'lure distractor exists vs. none' and 'transition from boundary node to within node vs. from boundary to boundary'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RTsubj    = zeros(nTrials, subLen);
RTtransTp = zeros(subLen, 2);
RTtransTp_type = zeros(subLen, 2, 2); % the 1st and 2nd '2' denote 'within vs. between cluster transition' and 'random and hamiltonian' walk
RTtransTp_Blc  = zeros(subLen, 2, nBlock);
RTtransTp_type_Blc = cell(1, 2); % random vs. hamiltonian walks
RTtransTp_lure = zeros(subLen, 2, 2); % the 1st and 2nd '2' denote 'lure distractor exists vs. none' and 'transition from boundary node to within node vs. from boundary to boundary'
lureTrialLen   = zeros(subLen, 2 * 2);
RTtransTp_type_lure = zeros(subLen, 2, 2, 2); % the last 2 denotes 'random and hamiltonian' walk
lme_subj = [];
for SubIdx = 1 : subLen
    subjBv      = subj_listBv{SubIdx};
    subTit      = subjLab{SubIdx};
    subjDir     = [bhvDataDir, subjBv, '/'];
    if isequal(ExpWord, 'cueEyeTrack')
        load([subjDir, subjBv, 'clusterResult_BlocksMouse.mat'], 'clusterResult');
    else
        load([subjDir, subjBv, 'clusterResult_Blocks.mat'], 'clusterResult');
    end
    trials_Col   = (1 : 1 : nTrials)';
    respRT_Col   = clusterResult(:, respRT);
    errorIdx_Col = clusterResult(:, errorIdx);
    timeOIdx_Col = clusterResult(:, timeOIdx);
    bonus_Col    = clusterResult(:, bonusCn);
    rndHam_Col   = clusterResult(:, trialTp); 
    fast_Col     = clusterResult(:, fastIdx); 
    corrIdx      = find(errorIdx_Col ~= 1);
    respRT_tmp   = respRT_Col;
    respRT_tmp(errorIdx_Col == 1) = nan;
    RTsubj(:, SubIdx) = respRT_tmp;
    %%% ratio of correct responses through predictions acorss 7 blocks
    for iBlock = 1 : nBlock
        ttBlcIdx  = (trialsInBlc(iBlock, 1) : trialsInBlc(iBlock, 2))';
        trialsBlc = length(ttBlcIdx);
        if isequal(ExpWord, 'cue') || isequal(ExpWord, 'cueEyeTrack') || isequal(ExpWord, 'cueRand')
            predRatio(SubIdx, iBlock) = length(find(respRT_Col(ttBlcIdx) <= grayoffT & errorIdx_Col(ttBlcIdx) ~= 1 & fast_Col(ttBlcIdx) ~= 1)) / trialsBlc; %% equaling to the ratio of 100 bonus
        elseif isequal(ExpWord, 'explicit')
            predRatio(SubIdx, iBlock) = length(find(bonus_Col(ttBlcIdx) == marioWin & errorIdx_Col(ttBlcIdx) ~= 1 & fast_Col(ttBlcIdx) ~= 1)) / trialsBlc; 
        end
        errRatio(SubIdx, iBlock)  = length(find(errorIdx_Col(ttBlcIdx) == 1 | timeOIdx_Col(ttBlcIdx) == 1)) / trialsBlc;
        fastRatio(SubIdx, iBlock) = length(find(fast_Col(ttBlcIdx) == 1)) / trialsBlc;
    end
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
    dtNum = arrayfun(@(x) length(find(dt_nodes(x, :)~=0)), (1 : nTrials)');
    
    %% baseline accuracy for each subject
    baseAcc(SubIdx, 1) = 1/2 * ((length(find(dtNum == 2))) / length(dtNum)) + ...
                         1/3 * ((length(find(dtNum == 3))) / length(dtNum)) + ...
                         1/4 * ((length(find(dtNum == 4))) / length(dtNum));
    
    %% accuracy in within-cluster and between-cluster transitions
    % added by rxj @ 08/30/2021
    transTpAcc(SubIdx, 1) = length(find(transStyle == 1 & errorIdx_Col ~= 1 & fast_Col ~= 1)) / length(find(transStyle == 1)); 
    transTpAcc(SubIdx, 2) = length(find(transStyle ~= 1 & errorIdx_Col ~= 1 & fast_Col ~= 1)) / length(find(transStyle ~= 1)); 
    
    %% Accuracy comparison: lure distractor exists vs. none
    % added by rxj @ 08/30/2021
    transTpAcc_lure(SubIdx, 1, 1) = length(find(nodesLabel == 1 & transStyle == 1 & lureIn == 1 & errorIdx_Col ~= 1 & fast_Col ~= 1)) / length(find(nodesLabel == 1 & transStyle == 1 & lureIn == 1));
    transTpAcc_lure(SubIdx, 1, 2) = length(find(nodesLabel == 1 & transStyle ~= 1 & lureIn == 1 & errorIdx_Col ~= 1 & fast_Col ~= 1)) / length(find(nodesLabel == 1 & transStyle ~= 1 & lureIn == 1));
    transTpAcc_lure(SubIdx, 2, 1) = length(find(nodesLabel == 1 & transStyle == 1 & lureIn == 0 & errorIdx_Col ~= 1 & fast_Col ~= 1)) / length(find(nodesLabel == 1 & transStyle == 1 & lureIn == 0));
    transTpAcc_lure(SubIdx, 2, 2) = length(find(nodesLabel == 1 & transStyle ~= 1 & lureIn == 0 & errorIdx_Col ~= 1 & fast_Col ~= 1)) / length(find(nodesLabel == 1 & transStyle ~= 1 & lureIn == 0));
    
    %% Comparison: lure distractor exists vs. none
    %RTtransTp_lure = zeros(subLen, 2, 2); % the 1st and 2nd '2' denote 'lure distractor exists vs. none' and 'transition from boundary node to within node vs. from boundary to boundary'
    RTtransTp_lure(SubIdx, 1, 1) = mean(respRT_Col(find(nodesLabel == 1 & transStyle == 1 & lureIn == 1 & errorIdx_Col ~= 1 & fast_Col ~= 1)));
    RTtransTp_lure(SubIdx, 1, 2) = mean(respRT_Col(find(nodesLabel == 1 & transStyle ~= 1 & lureIn == 1 & errorIdx_Col ~= 1 & fast_Col ~= 1)));
    RTtransTp_lure(SubIdx, 2, 1) = mean(respRT_Col(find(nodesLabel == 1 & transStyle == 1 & lureIn == 0 & errorIdx_Col ~= 1 & fast_Col ~= 1)));
    RTtransTp_lure(SubIdx, 2, 2) = mean(respRT_Col(find(nodesLabel == 1 & transStyle ~= 1 & lureIn == 0 & errorIdx_Col ~= 1 & fast_Col ~= 1)));
    
    lureTrialLen(SubIdx, 1) = length(find(nodesLabel == 1 & transStyle == 1 & lureIn == 1 & errorIdx_Col ~= 1 & fast_Col ~= 1));
    lureTrialLen(SubIdx, 2) = length(find(nodesLabel == 1 & transStyle ~= 1 & lureIn == 1 & errorIdx_Col ~= 1 & fast_Col ~= 1));
    lureTrialLen(SubIdx, 3) = length(find(nodesLabel == 1 & transStyle == 1 & lureIn == 0 & errorIdx_Col ~= 1 & fast_Col ~= 1));
    lureTrialLen(SubIdx, 4) = length(find(nodesLabel == 1 & transStyle ~= 1 & lureIn == 0 & errorIdx_Col ~= 1 & fast_Col ~= 1));
    
    %% Comparison: lure distractor exists vs. none, in random and hamiltonian path
    %RTtransTp_type_lure = zeros(subLen, 2, 2, 2); % the last 2 denotes 'random and hamiltonian' walk
    %%% wtihin vs. between cluster transition RTs separately for random and
    %%% hamiltonian walk
    RTtransTp_type_lure(SubIdx, 1, 1, 1) = mean(respRT_Col(find(nodesLabel == 1 & transStyle == 1 & lureIn == 1 & errorIdx_Col ~= 1 & fast_Col ~= 1 & rndHam_Col == 1))); %% rndHam_Col == 1: random trial
    RTtransTp_type_lure(SubIdx, 1, 2, 1) = mean(respRT_Col(find(nodesLabel == 1 & transStyle ~= 1 & lureIn == 1 & errorIdx_Col ~= 1 & fast_Col ~= 1 & rndHam_Col == 1))); 
    RTtransTp_type_lure(SubIdx, 1, 1, 2) = mean(respRT_Col(find(nodesLabel == 1 & transStyle == 1 & lureIn == 1 & errorIdx_Col ~= 1 & fast_Col ~= 1 & rndHam_Col == 2))); %% rndHam_Col == 2: hamiltonian trial
    RTtransTp_type_lure(SubIdx, 1, 2, 2) = mean(respRT_Col(find(nodesLabel == 1 & transStyle ~= 1 & lureIn == 1 & errorIdx_Col ~= 1 & fast_Col ~= 1 & rndHam_Col == 2))); 
    
    RTtransTp_type_lure(SubIdx, 2, 1, 1) = mean(respRT_Col(find(nodesLabel == 1 & transStyle == 1 & lureIn == 0 & errorIdx_Col ~= 1 & fast_Col ~= 1 & rndHam_Col == 1)));
    RTtransTp_type_lure(SubIdx, 2, 2, 1) = mean(respRT_Col(find(nodesLabel == 1 & transStyle ~= 1 & lureIn == 0 & errorIdx_Col ~= 1 & fast_Col ~= 1 & rndHam_Col == 1)));
    RTtransTp_type_lure(SubIdx, 2, 1, 2) = mean(respRT_Col(find(nodesLabel == 1 & transStyle == 1 & lureIn == 0 & errorIdx_Col ~= 1 & fast_Col ~= 1 & rndHam_Col == 2)));
    RTtransTp_type_lure(SubIdx, 2, 2, 2) = mean(respRT_Col(find(nodesLabel == 1 & transStyle ~= 1 & lureIn == 0 & errorIdx_Col ~= 1 & fast_Col ~= 1 & rndHam_Col == 2))); 
    
    %%
    RTtransTp(SubIdx, 1) = mean(respRT_Col(find(transStyle == 1 & errorIdx_Col ~= 1 & fast_Col ~= 1))); %% Pay attention: errorIdx_Col includes error and timeout trials 
    RTtransTp(SubIdx, 2) = mean(respRT_Col(find(transStyle ~= 1 & errorIdx_Col ~= 1 & fast_Col ~= 1)));
    %%% wtihin vs. between cluster transition RTs separately for random and
    %%% hamiltonian walk
    RTtransTp_type(SubIdx, 1, 1) =  mean(respRT_Col(find(transStyle == 1 & errorIdx_Col ~= 1 & fast_Col ~= 1 & rndHam_Col == 1))); %% rndHam_Col == 1: random trial
    RTtransTp_type(SubIdx, 1, 2) =  mean(respRT_Col(find(transStyle == 1 & errorIdx_Col ~= 1 & fast_Col ~= 1 & rndHam_Col == 2))); %% rndHam_Col == 2: hamiltonian trial
    RTtransTp_type(SubIdx, 2, 1) =  mean(respRT_Col(find(transStyle ~= 1 & errorIdx_Col ~= 1 & fast_Col ~= 1 & rndHam_Col == 1)));
    RTtransTp_type(SubIdx, 2, 2) =  mean(respRT_Col(find(transStyle ~= 1 & errorIdx_Col ~= 1 & fast_Col ~= 1 & rndHam_Col == 2)));
    %%% within vs. between cluster transition RTs across blocks
    for iBlock = 1 : nBlock
        ttBlcIdx = (trialsInBlc(iBlock, 1) : trialsInBlc(iBlock, 2))';
        RTtransTp_Blc(SubIdx, 1, iBlock) = mean(respRT_Col(find((trials_Col >= ttBlcIdx(1) & trials_Col <= ttBlcIdx(end)) & transStyle == 1 & errorIdx_Col ~= 1 & fast_Col ~= 1)));
        RTtransTp_Blc(SubIdx, 2, iBlock) = mean(respRT_Col(find((trials_Col >= ttBlcIdx(1) & trials_Col <= ttBlcIdx(end)) & transStyle ~= 1 & errorIdx_Col ~= 1 & fast_Col ~= 1)));
    end
    %%% within vs. between cluster transition RTs across blocks seperately for random and hamiltonian walk
    for rOh = 1 : 2 % RTtransTp_type_Blc = cell(1, 2);
        for iBlock = 1 : nBlock
            ttBlcIdx = (trialsInBlc(iBlock, 1) : trialsInBlc(iBlock, 2))';
            RTtransTp_type_Blc{rOh}(SubIdx, 1, iBlock) = mean(respRT_Col(find((trials_Col >= ttBlcIdx(1) & trials_Col <= ttBlcIdx(end)) & transStyle == 1 & errorIdx_Col ~= 1 & fast_Col ~= 1 & rndHam_Col == rOh)));
            RTtransTp_type_Blc{rOh}(SubIdx, 2, iBlock) = mean(respRT_Col(find((trials_Col >= ttBlcIdx(1) & trials_Col <= ttBlcIdx(end)) & transStyle ~= 1 & errorIdx_Col ~= 1 & fast_Col ~= 1 & rndHam_Col == rOh)));
        end
    end
    %%% table for the lme test: RT, 
    lme_i = [];
    recency_nodes = zeros(nTrials, 1);
    for iT = 1 : nTrials
        rec_iT = find(from_nodes(iT) == from_nodes(1 : (iT - 1)));
        if isempty(rec_iT)
            recency_nodes(iT) = nodeNum;
        else 
            recency_nodes(iT) = iT - rec_iT(end);
        end
    end
    is_BdNode = zeros(nTrials, 1); %% whether boundary nodes
    is_BdNode(find(from_nodes == 1 | from_nodes == 5 | from_nodes == 6 | from_nodes == 10 | from_nodes == 11 | from_nodes == 15)) = 1;
    lme_i = [repmat(SubIdx, nTrials, 1), trials_Col, respRT_Col, from_nodes, to_nodes, dtNum, is_BdNode, recency_nodes, transStyle, rndHam_Col];
    lme_i = lme_i(corrIdx, :);
    lme_i(:, 2) = log(lme_i(:, 2)); %% log(trial)
    lme_subj = [lme_subj; lme_i];
end

%% color settings
colorSet = [249, 183, 176; ...
            138, 170, 51; ...
            84, 185, 211; ...
            248, 218, 172; ...
            184, 204, 225; ...
            210, 234, 200; ...
            198, 127, 192; ...
            219, 204, 226] ./ [255, 255, 255];
colorSubj = hsv(subLen); 

%% Learning effect: ratio of error or timeout trials / too fast ratio
[erAvg, erSem] = Mean_and_Se(fastRatio, 1);
%%% (1) plot the predicted response (or full bonus=100) across the 7 blocks
LineSty = '-';
% for SubIdx = 1 : subLen
%     subTit = subjLab{SubIdx};
%     figure('Position', [100 100 400 300]), clf;
%     plot(1 : 1 : nBlock, fastRatio(SubIdx, :), 'Marker', '.', 'MarkerSize', 40, 'Color', [0.4, 0.4, 0.4], 'LineStyle', LineSty, 'LineWidth', 4); hold on;
%     set(gca, 'FontSize', 20, 'FontWeight', 'Bold', 'LineWidth', 2);
%     set(gca, 'XTick', 1 : 1 : nBlock, 'XTickLabel', '');
%     xlim([1, nBlock]);
%     box off;
%     title(subTit, 'FontSize', 20);
% end
figure('Position', [100 100 400 300]), clf;
errorbar(1 : 1 : nBlock, erAvg, erSem, 'Color', 'k', 'LineStyle', '-', 'LineWidth', 4); hold on;
%plot(1 : 1 : nBlock, prAvg, 'Marker', '.', 'MarkerSize', 40, 'Color', [0, 0, 0], 'LineStyle', '-', 'LineWidth', 4); hold on;
set(gca, 'FontSize', 20, 'FontWeight', 'Bold', 'LineWidth', 2);
set(gca, 'XTick', 1 : 1 : nBlock, 'XTickLabel', 1 : 1 : nBlock);
xlim([1, nBlock]);
box off;
xlabel('Blocks', 'FontSize', 20);
ylabel('Fast ratio', 'FontSize', 20);


%% Learning effect (individual participant): the correct response time decreases with the trials
nBins = nTrials / HamInter;
for SubIdx = 1 : subLen
    subTit = subjLab{SubIdx};
    RT_i   = RTsubj(:, SubIdx);
    RT_bins = zeros(nBins, 1);
    for iB = 1 : nBins
        iB_id = (iB - 1) * HamInter + 1 : iB * HamInter;
        RT_bins(iB) = nanmean(RT_i(iB_id));
    end
    figure('Position', [100 100 400 300]), clf;
    %plot(1 : length(RT_i), RT_i, 'Marker', '.', 'MarkerSize', 10, 'Color', [0.4, 0.4, 0.4], 'LineStyle', 'none'); hold on;
    plot(1 : 1 : nBins, RT_bins, 'Marker', '.', 'MarkerSize', 20, 'Color', [0.4, 0.4, 0.4], 'LineStyle', '-', 'LineWidth', 3); hold on;
    set(gca, 'FontSize', 20, 'FontWeight', 'Bold', 'LineWidth', 2);
    set(gca, 'XTick', 0 : 10 : nBins, 'XTickLabel', '');
    ylim([0.5, 2.5]);
    box off;
    xlabel('Trials', 'FontSize', 20);
    ylabel('Correct RTs (s)', 'FontSize', 20);
    title(subTit, 'FontSize', 20);
end

%% Learning effect: the predicted response ratio increases with the blocks
[prAvg, prSem] = Mean_and_Se(predRatio, 1);
%%% (1) plot the predicted response (or full bonus=100) across the 7 blocks
figure('Position', [100 100 1500 800]), clf;
hold on;
if isequal(ExpWord, 'cue')
    ha = tight_subplot(4, 6, [.05 .05], [.1 .05], [.05 .05]); %% [ha, pos] = tight_subplot(Nh, Nw, gap, marg_h, marg_w)
elseif isequal(ExpWord, 'explicit')
    ha = tight_subplot(4, 6, [.05 .05], [.1 .05], [.05 .05]); %% [ha, pos] = tight_subplot(Nh, Nw, gap, marg_h, marg_w)
elseif isequal(ExpWord, 'cueRand')
    ha = tight_subplot(4, 6, [.05 .05], [.1 .05], [.05 .05]); 
elseif isequal(ExpWord, 'cueEyeTrack')
    ha = tight_subplot(3, 6, [.05 .05], [.1 .05], [.05 .05]);
end

ha_i = 1;
LineSty = '-';
for SubIdx = 1 : subLen
    subTit = subjLab{SubIdx};
    %figure('Position', [100 100 400 300]), clf;
    axes(ha(ha_i));
    plot(1 : 1 : nBlock, predRatio(SubIdx, :), 'Marker', '.', 'MarkerSize', 40, 'Color', [0.4, 0.4, 0.4], 'LineStyle', LineSty, 'LineWidth', 4); hold on;
    set(gca, 'FontSize', 16, 'FontWeight', 'Bold', 'LineWidth', 2);
    set(gca, 'XTick', 1 : 1 : nBlock, 'XTickLabel', '');
    xlim([1, nBlock]);
    box off;
    title(subTit, 'FontSize', 16);
    axis('square');
    ha_i = ha_i + 1;
end
%%% *************plotting the average across subjects
figure('Position', [100 100 400 300]), clf;
errorbar(1 : 1 : nBlock, prAvg, prSem, 'Color', 'k', 'LineStyle', '-', 'LineWidth', 4); hold on;
%plot(1 : 1 : nBlock, prAvg, 'Marker', '.', 'MarkerSize', 40, 'Color', [0, 0, 0], 'LineStyle', '-', 'LineWidth', 4); hold on;
set(gca, 'FontSize', 20, 'FontWeight', 'Bold', 'LineWidth', 2);
set(gca, 'XTick', 1 : 1 : nBlock, 'XTickLabel', 1 : 1 : nBlock);
xlim([1, nBlock]);
box off;
axis('square');
xlabel('Blocks', 'FontSize', 20);
ylabel('Prediction ratio', 'FontSize', 20);
%%
figure('Position', [100 100 400 300]), clf;
errorbar(1 : 1 : nBlock, prAvg, prSem, 'Color', 'k', 'LineStyle', '-', 'LineWidth', 4); hold on;
%plot(1 : 1 : nBlock, prAvg, 'Marker', '.', 'MarkerSize', 40, 'Color', [0, 0, 0], 'LineStyle', '-', 'LineWidth', 4); hold on;
set(gca, 'FontSize', 20, 'FontWeight', 'Bold', 'LineWidth', 2);
set(gca, 'XTick', 1 : 1 : nBlock, 'XTickLabel', 1 : 1 : nBlock);
xlim([1, nBlock]);
box off;
%axis('square');
% xlabel('Blocks', 'FontSize', 20);
% ylabel('Prediction ratio', 'FontSize', 20);



%% Learning effect: cross-cluster RT vs. within-cluster RT
% (1) random + hamiltonian
[RTavg_all, RTsem_all] = Mean_and_Se(RTtransTp, 1);
figure('Position', [100 100 300 200]), clf;
for iB = 1 : 2
    b = bar(iB, RTavg_all(iB), 0.45, 'LineStyle', '-', 'LineWidth', 2); hold on;
    b.FaceColor = colorSet(iB, :);
    errorbar(iB, RTavg_all(iB), RTsem_all(iB), 'Color', 'k', 'LineStyle', 'none', 'LineWidth', 2); hold on;
end
set(gca, 'FontSize', 20, 'FontWeight', 'Bold', 'LineWidth', 2);
set(gca, 'XTick', '', 'XTickLabel', '');
xlim([1-0.5, 2+0.5]);
box off;
[h, p, ci, stats] = ttest(RTtransTp(:, 2), RTtransTp(:, 1))
%%% the difference
figure('Position', [100 100 150 200]), clf;
[RTavg_allDif, RTsem_allDif] = Mean_and_Se(RTtransTp(:, 2) - RTtransTp(:, 1), 1);
b = bar(1, RTavg_allDif, 0.45, 'LineStyle', '-', 'LineWidth', 2); hold on;
b.FaceColor = colorSet(iB + 1, :);
errorbar(1, RTavg_allDif, RTsem_allDif, 'Color', 'k', 'LineStyle', 'none', 'LineWidth', 2); hold on;
set(gca, 'FontSize', 20, 'FontWeight', 'Bold', 'LineWidth', 2);
set(gca, 'XTick', '', 'XTickLabel', '');
xlim([1-0.5, 1+0.5]);
box off;

%% (2) random & (3) hamiltonian
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
iExp = 1;
colorTmp = colorExp(iExp, :);

[RTavg_rOh, RTsem_rOh] = Mean_and_Se(RTtransTp_type, 1);
RTavg_rOh = squeeze(RTavg_rOh);
RTsem_rOh = squeeze(RTsem_rOh);% RTtransTp_type = zeros(subLen, 2, 2); % the 1st and 2nd '2' denote 'within vs. between cluster transition' and 'random and hamiltonian' walk
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
    plot(barPos_i, RTtransTp_type(:, :, iB), '.-', 'LineWidth', 2, 'MarkerSize', 10, 'Color', [0.6, 0.6, 0.6]); hold on; % individuals
    plot(barPos_i, RTavg_rOh(:, iB), 'Marker', '.', 'MarkerSize', 10, 'Color', colorTmp, 'LineStyle', 'none'); hold on;
    errorbar(barPos_i(1), RTavg_rOh(1, iB), RTsem_rOh(1, iB), 'Color', colorTmp, 'Marker', '.', 'LineStyle', '-', 'LineWidth', 4, 'MarkerSize', 10);
    errorbar(barPos_i(2), RTavg_rOh(2, iB), RTsem_rOh(2, iB), 'Color', colorTmp, 'Marker', '.', 'LineStyle', '-', 'LineWidth', 4, 'MarkerSize', 10);
%     for i_rOh = 1 : 2
%         %%% bar plot
%         b = bar(barPos_i(i_rOh), RTavg_rOh(i_rOh, iB), 0.45, 'LineStyle', LineBar, 'LineWidth', 2); hold on;
%         b.FaceColor = colorSet(i_rOh, :);
%         errorbar(barPos_i(i_rOh), RTavg_rOh(i_rOh, iB), RTsem_rOh(i_rOh, iB), 'Color', 'k', 'LineStyle', 'none', 'LineWidth', 2); hold on;
%     end
end
set(gca, 'FontSize', 16, 'FontWeight', 'Bold', 'LineWidth', 2);
set(gca, 'XTick', '', 'XTickLabel', '');
xlim([0.9, barPos(end, end)+0.1]);
ylim([0.5, 2]);
box off;

%% the difference
[RTavg_rOhDif, RTsem_rOhDif] = Mean_and_Se(RTtransTp_type(:, 2, :) - RTtransTp_type(:, 1, :), 1);
RTavg_rOhDif = squeeze(RTavg_rOhDif);
RTsem_rOhDif = squeeze(RTsem_rOhDif);
figure('Position', [100 100 300 200]), clf;
for iB = 1 : 2
    if iB == 1     %% random walk
        LineBar = '-';
    elseif iB == 2 %% hamiltonian walk
        LineBar = ':';
    end
    b = bar(iB, RTavg_rOhDif(iB), 0.45, 'LineStyle', LineBar, 'LineWidth', 2); hold on;
    b.FaceColor = colorSet(3, :);
    errorbar(iB, RTavg_rOhDif(iB), RTsem_rOhDif(iB), 'Color', 'k', 'LineStyle', 'none', 'LineWidth', 2); hold on;
end
set(gca, 'FontSize', 20, 'FontWeight', 'Bold', 'LineWidth', 2);
set(gca, 'XTick', '', 'XTickLabel', '');
xlim([1-0.5, 2+0.5]);
box off;
[h_r, p_r, ci_r, stats_r] = ttest(RTtransTp_type(:, 2, 1), RTtransTp_type(:, 1, 1))
[h_h, p_h, ci_h, stats_h] = ttest(RTtransTp_type(:, 2, 2), RTtransTp_type(:, 1, 2))
%%% whether there is significant difference between random and hamiltonian
%%% walk
[h_rOh, p_rOh, ci_rOh, stats_rOh] = ttest(RTtransTp_type(:, 2, 2)-RTtransTp_type(:, 1, 2), RTtransTp_type(:, 2, 1)-RTtransTp_type(:, 1, 1))

%% anova on RTs, with variables of random or hamiltonian, and within- or between-cluster transitions
% RTtransTp_type = zeros(subLen, 2, 2);  % the 1st and 2nd '2' denote 'within vs. between cluster transition' and 'random and hamiltonian' walk
%%% repeated measurement ANOVA for RTs
%  2 (random vs. hamiltonian) by 2 (within- vs. between-cluster transitions) repeated-measures ANOVAs
nWalk  = 2;
nTran  = 2;
ID     = reshape(repmat((1 : 1 : subLen), nWalk * nTran, 1), subLen * nWalk * nTran, 1);
walkCol  = repmat((1 : 1 : nWalk)', subLen, nTran);
walkCol  = reshape(walkCol', size(walkCol, 1)*size(walkCol, 2), 1);
transCol = repmat([1; 2], subLen * nWalk, 1); 
model = [1, 0, 0;... % main effect: subject
         0, 1, 0;... % main effect: walks
         0, 0, 1; ...% main effect: transitions
         0, 1, 1]; 
RTtest = reshape(RTtransTp_type, [subLen, nTran * nWalk]);
RTtest = reshape(RTtest', [subLen * nTran * nWalk, 1]);
[p_RTs, tbl_RTs, stats_RTs] = anovan(RTtest, {ID, walkCol, transCol}, 'model', model, 'random', [1], 'sstype', 3, 'varnames', {'Subjects', 'WalksType', 'TransType'}); % 'model' as 'interaction' is completely the same as model1


%% Accuracy in random vs. hamiltonian walks
[accRatio_rOh, accRatio_tp] = ClusterRep_v2_predRatio(bhvDataDir, subj_listBv, subLen, errorIdx, fastIdx, trialTp, objTgt, objResp, nTrials, clsDef);
%%% Figure 1&2: only trial types (random or hamiltonian) or transition types (within or between)
%%% Figure 3: combinations of trial types & transition types
for i = 1 : 2
    if i == 1
        accMat   = accRatio_rOh;
        titWor   = 'Random or Hamiltonian: ';
        colorTmp = [colorSet(4, :); colorSet(4, :)];
        LineSty  = {'-', ':'};
    elseif i == 2
        accMat   = accRatio_tp;
        titWor   = 'Within- or between-transition: ';
        colorTmp = colorSet(1 : 2, :);
        LineSty  = {'-', '-'};
    end
    [accAvg, accSem] = Mean_and_Se(accMat, 1); % accRatio_rOh = zeros(subLen, 2); %% correct response in random and hamiltonian walks
    figure('Position', [100 100 300 200]), clf;
    for iB = 1 : size(accAvg, 2)
        b = bar(iB, accAvg(iB), 0.45, 'LineStyle', LineSty{iB}, 'LineWidth', 2); hold on;
        b.FaceColor = colorTmp(iB, :);
        errorbar(iB, accAvg(iB), accSem(iB), 'Color', 'k', 'LineStyle', 'none', 'LineWidth', 2); hold on;
    end
    figureSet('', '', 0.5, size(accMat, 2) + 0.5);
    disp(titWor);
    [h, p, ci, stats] = ttest(accMat(:, 2), accMat(:, 1))

end

%% Accuracy in within- & between- transition, also random & hamiltonian walks
[~, ~, accRatio_tp_rOh] = ClusterRep_v2_predRatio(bhvDataDir, subj_listBv, subLen, errorIdx, fastIdx, trialTp, objTgt, objResp, nTrials, clsDef);
accMat   = accRatio_tp_rOh; %% accRatio_tp_rOh = zeros(subLen, 2, 2); %% the 1st and 2nd '2' denote ' within-cluster vs. between-cluster transition' and 'random & hamiltonian walks'
titWor   = '2 walks & 2 transitions: ';
colorTmp = colorSet(1 : 2, :);
LineSty  = {'-', ':'; '-', ':'};
[accAvg, accSem] = Mean_and_Se(accMat, 1);
accAvg = squeeze(accAvg);
accSem = squeeze(accSem);

figure('Position', [100 100 400 200]), clf;
barPos = [1, 1.5; 2.5, 3];
for iB = 1 : 2
    barPos_i = barPos(iB, :);
    if iB == 1     %% random walk
        LineBar = '-';
    elseif iB == 2 %% hamiltonian walk
        LineBar = ':';
    end
    for i_rOh = 1 : 2
        b = bar(barPos_i(i_rOh), accAvg(i_rOh, iB), 0.45, 'LineStyle', LineBar, 'LineWidth', 2); hold on;
        b.FaceColor = colorTmp(i_rOh, :);
        errorbar(barPos_i(i_rOh), accAvg(i_rOh, iB), accSem(i_rOh, iB), 'Color', 'k', 'LineStyle', 'none', 'LineWidth', 2); hold on;
    end
end
set(gca, 'FontSize', 20, 'FontWeight', 'Bold', 'LineWidth', 2);
set(gca, 'XTick', '', 'XTickLabel', '');
xlim([0.5, barPos(end, end)+0.5]);
box off;
%%% the difference between random and hamiltonian walks separately for
%%% within- and between-cluster transitions
[accAvg_Dif, accSem_Dif] = Mean_and_Se(accRatio_tp_rOh(:, 2, :) - accRatio_tp_rOh(:, 1, :), 1);
accAvg_Dif = squeeze(accAvg_Dif);
accSem_Dif = squeeze(accSem_Dif);
figure('Position', [100 100 300 200]), clf;
for iB = 1 : 2
    if iB == 1     %% random walk
        LineBar = '-';
    elseif iB == 2 %% hamiltonian walk
        LineBar = ':';
    end
    b = bar(iB, accAvg_Dif(iB), 0.45, 'LineStyle', LineBar, 'LineWidth', 2); hold on;
    b.FaceColor = colorSet(3, :);
    errorbar(iB, accAvg_Dif(iB), accSem_Dif(iB), 'Color', 'k', 'LineStyle', 'none', 'LineWidth', 2); hold on;
end
set(gca, 'FontSize', 20, 'FontWeight', 'Bold', 'LineWidth', 2);
set(gca, 'XTick', '', 'XTickLabel', '');
xlim([1-0.5, 2+0.5]);
box off;
[h_r, p_r, ci_r, stats_r] = ttest(accRatio_tp_rOh(:, 2, 1), accRatio_tp_rOh(:, 1, 1)) %% between- vs. within-cluster transitions in random walks
[h_h, p_h, ci_h, stats_h] = ttest(accRatio_tp_rOh(:, 2, 2), accRatio_tp_rOh(:, 1, 2)) %% between- vs. within-cluster transitions in hamiltonian walks
%%% whether there is significant difference between random and hamiltonian
%%% walk
[h_rOh, p_rOh, ci_rOh, stats_rOh] = ttest(accRatio_tp_rOh(:, 2, 2)-accRatio_tp_rOh(:, 1, 2), accRatio_tp_rOh(:, 2, 1)-accRatio_tp_rOh(:, 1, 1))


%% anova on correct response ratios, with variables of random or hamiltonian, and within- or between-cluster transitions
% accMat   = accRatio_tp_rOh; %% accRatio_tp_rOh = zeros(subLen, 2, 2); %% the 1st and 2nd '2' denote ' within-cluster vs. between-cluster transition' and 'random & hamiltonian walks'
%%% repeated measurement ANOVA for corect response ratios
%  2 (random vs. hamiltonian) by 2 (within- vs. between-cluster transitions) repeated-measures ANOVAs
nWalk  = 2;
nTran  = 2;
ID     = reshape(repmat((1 : 1 : subLen), nWalk * nTran, 1), subLen * nWalk * nTran, 1);
walkCol  = repmat((1 : 1 : nWalk)', subLen, nTran);
walkCol  = reshape(walkCol', size(walkCol, 1)*size(walkCol, 2), 1);
transCol = repmat([1; 2], subLen * nWalk, 1); 
model = [1, 0, 0;... % main effect: subject
         0, 1, 0;... % main effect: walks
         0, 0, 1; ...% main effect: transitions
         0, 1, 1]; 
RTtest = reshape(accRatio_tp_rOh, [subLen, nTran * nWalk]);
RTtest = reshape(RTtest', [subLen * nTran * nWalk, 1]);
[p_ratio, tbl_ratio, stats_ratio] = anovan(RTtest, {ID, walkCol, transCol}, 'model', model, 'random', [1], 'sstype', 3, 'varnames', {'Subjects', 'WalksType', 'TransType'}); % 'model' as 'interaction' is completely the same as model1

%% correlation between accuracy and RTs across subjects
% added by rxj @ 04/22/2021
[accRatio_rOh, accRatio_tp, accRatio_tp_rOh] = ClusterRep_v2_predRatio(bhvDataDir, subj_listBv, subLen, errorIdx, fastIdx, trialTp, objTgt, objResp, nTrials, clsDef);
% RTtransTp = zeros(subLen, 2);
% RTtransTp_type = zeros(subLen, 2, 2);
[r_tpW, p_tpW] = corr(RTtransTp(:, 1), accRatio_tp(:, 1))
[r_tpB, p_tpB] = corr(RTtransTp(:, 2), accRatio_tp(:, 2))
[r_tpBW, p_tpBW] = corr(RTtransTp(:, 2)-RTtransTp(:, 1), accRatio_tp(:, 2)-accRatio_tp(:, 1))
figure('Position', [100 100 600 200]), clf;
subplot(1, 2, 1);
plot(RTtransTp(:, 1), accRatio_tp(:, 1), 'Marker', '.', 'MarkerSize', 18, 'Color', colorSet(1, :), 'MarkerFaceColor', colorSet(1, :), 'LineStyle', 'none'); hold on;
plot(RTtransTp(:, 2), accRatio_tp(:, 2), 'Marker', '.', 'MarkerSize', 18, 'Color', colorSet(2, :), 'MarkerFaceColor', colorSet(2, :), 'LineStyle', 'none'); hold on;
set(gca, 'FontSize', 20, 'FontWeight', 'Bold', 'LineWidth', 2);
%set(gca, 'XTick', '', 'XTickLabel', '');
box off;
subplot(1, 2, 2);
plot(RTtransTp(:, 2)-RTtransTp(:, 1), accRatio_tp(:, 2)-accRatio_tp(:, 1), 'Marker', '.', 'MarkerSize', 18, 'Color', colorSet(3, :), 'MarkerFaceColor', colorSet(3, :), 'LineStyle', 'none'); hold on;
set(gca, 'FontSize', 20, 'FontWeight', 'Bold', 'LineWidth', 2);
%set(gca, 'XTick', '', 'XTickLabel', '');
box off;


%% with lure distractors vs. without
%% which matrix
flg = 3;
if flg == 1
    RTmat = RTtransTp_lure; %% mixture of random and hamiltonian path
elseif flg == 2
    RTmat = RTtransTp_type_lure(:, :, :, 1); %% only random path
elseif flg == 3
    RTmat = RTtransTp_type_lure(:, :, :, 2); %% only hamiltonian path
end
%% with lure distractors vs. without
[RTavg_rOh, RTsem_rOh] = Mean_and_Se(RTmat, 1);
RTavg_rOh = squeeze(RTavg_rOh);
RTsem_rOh = squeeze(RTsem_rOh); % the 1st and 2nd '2' denote 'lure distractor exists vs. none' and 'transition from boundary node to within node vs. from boundary to boundary'
figure('Position', [100 100 400 200]), clf;
barPos = [1, 1.5; 2.5, 3];
for iB = 1 : 2
    barPos_i = barPos(iB, :);
    if iB == 1     %% random walk
        LineBar = '-';
    elseif iB == 2 %% hamiltonian walk
        LineBar = ':';
    end
    for i_rOh = 1 : 2
        b = bar(barPos_i(i_rOh), RTavg_rOh(i_rOh, iB), 0.45, 'LineStyle', LineBar, 'LineWidth', 2); hold on;
        b.FaceColor = colorSet(i_rOh, :);
        errorbar(barPos_i(i_rOh), RTavg_rOh(i_rOh, iB), RTsem_rOh(i_rOh, iB), 'Color', 'k', 'LineStyle', 'none', 'LineWidth', 2); hold on;
    end
end
set(gca, 'FontSize', 20, 'FontWeight', 'Bold', 'LineWidth', 2);
set(gca, 'XTick', '', 'XTickLabel', '');
xlim([0.5, barPos(end, end)+0.5]);
box off;
%%% the difference
[RTavg_rOhDif, RTsem_rOhDif] = Mean_and_Se(RTmat(:, 2, :) - RTmat(:, 1, :), 1);
RTavg_rOhDif = squeeze(RTavg_rOhDif);
RTsem_rOhDif = squeeze(RTsem_rOhDif);
figure('Position', [100 100 300 200]), clf;
for iB = 1 : 2
    if iB == 1     %% random walk
        LineBar = '-';
    elseif iB == 2 %% hamiltonian walk
        LineBar = ':';
    end
    b = bar(iB, RTavg_rOhDif(iB), 0.45, 'LineStyle', LineBar, 'LineWidth', 2); hold on;
    b.FaceColor = colorSet(3, :);
    errorbar(iB, RTavg_rOhDif(iB), RTsem_rOhDif(iB), 'Color', 'k', 'LineStyle', 'none', 'LineWidth', 2); hold on;
end
set(gca, 'FontSize', 20, 'FontWeight', 'Bold', 'LineWidth', 2);
set(gca, 'XTick', '', 'XTickLabel', '');
xlim([1-0.5, 2+0.5]);
box off;
[h_r, p_r, ci_r, stats_r] = ttest(RTmat(:, 2, 1), RTmat(:, 1, 1))
[h_h, p_h, ci_h, stats_h] = ttest(RTmat(:, 2, 2), RTmat(:, 1, 2))
%%% whether there is significant difference between random and hamiltonian
%%% walk
[h_rOh, p_rOh, ci_rOh, stats_rOh] = ttest(RTmat(:, 2, 2)-RTmat(:, 1, 2), RTmat(:, 2, 1)-RTmat(:, 1, 1))

%%% single subjects
figure('Position', [100 100 400 200]), clf;
for SubIdx = 1 : subLen
    for rOh = 1 : 2
        xLab = (rOh - 1) * 2 + 1 : rOh * 2;
        if rOh == 1
            LineSty = '-';
        elseif rOh == 2
            LineSty = ':';
        end
        plot(xLab, RTmat(SubIdx, :, rOh), 'Marker', '.', 'MarkerSize', 40, 'Color', colorSet(rOh + 3, :), 'LineStyle', LineSty, 'LineWidth', 2); hold on;
    end
end
set(gca, 'FontSize', 16, 'FontWeight', 'Bold', 'LineWidth', 2);
set(gca, 'XTick', 1 : 1 : 4, 'XTickLabel', '');
xlim([1-0.5, 4+0.5]);
box off;

%% Concatenating the respYes_trials across blocks for the subsequent use in Jupyter Lab
% write by rxj @ 03/16/2021
for SubIdx = 1 : subLen
    subjBv      = subj_listBv{SubIdx};
    subTit      = subjLab{SubIdx};
    subjDir     = [bhvDataDir, subjBv, '/'];
    %%% ratio of correct responses through predictions acorss 7 blocks
    mouseTraj_trials_blc = [];
    respYes_trials_blc   = [];
    for iBlock = 1 : nBlock
        load([subjDir, subjBv, 'mouseTraj_trials_in_Blc', num2str(iBlock), '.mat'], 'mouseTraj_trials');
        load([subjDir, subjBv, 'respYes_trials_in_Blc', num2str(iBlock), '.mat'], 'respYes_trials');
        mouseTraj_trials_blc = [mouseTraj_trials_blc; mouseTraj_trials];
        respYes_trials_blc   = [respYes_trials_blc; respYes_trials];
    end
    save([subjDir, subjBv, 'mouseTraj_trials_blc.mat'], 'mouseTraj_trials_blc');
    save([subjDir, subjBv, 'respYes_trials_blc.mat'], 'respYes_trials_blc');
end

%% Calculating the recency
% write by rxj @ 04/13/2021
recency_nodes_sub = zeros(nTrials, 2, subLen);
for SubIdx = 1 : subLen
    subjBv      = subj_listBv{SubIdx};
    subTit      = subjLab{SubIdx};
    subjDir     = [bhvDataDir, subjBv, '/'];
    load([subjDir, subjBv, 'clusterResult_Blocks.mat'], 'clusterResult');
    from_nodes = clusterResult(:, objTgt);

    recency_nodes = zeros(nTrials, 2);
    for iT = 1 : nTrials
        rec_iT = find(from_nodes(iT) == from_nodes(1 : (iT - 1)));
        if isempty(rec_iT)
            recency_nodes(iT, 1) = nodeNum;
            recency_nodes(iT, 2) = nodeNum;
            
        elseif ~isempty(rec_iT) && length(rec_iT) == 1
            recency_nodes(iT, 1) = iT - rec_iT(end);
            recency_nodes(iT, 2) = nodeNum;
            
        elseif ~isempty(rec_iT) && length(rec_iT) >= 2
            recency_nodes(iT, 1) = iT - rec_iT(end);
            recency_nodes(iT, 2) = iT - rec_iT(end - 1);
        end
    end
    recency_nodes_sub(:, :, SubIdx) = recency_nodes;
end








