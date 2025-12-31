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
ExpWord_List = {'ImplicitExp', 'ExplicitExp', 'ImplicitRandExp'};
ExpIdx       = 2;
ExpWord      = ExpWord_List{ExpIdx};
if isequal(ExpWord, 'ImplicitExp')
    subj_list = {'wsn_1_f_18', 'dy_2_f_22', 'haq_3_f_24', 'hry_4_f_20', 'zjx_5_m_20', 'yyq_6_f_18', 'zkw_7_m_18', 'zy_8_f_20', 'hys_9_m_20', 'cjj_10_m_18', ...
                 'dwq_11_f_22', 'ljl_12_m_20', 'jyx_13_m_19', 'zk_14_f_21', 'lsy_15_m_19', 'cjl_16_m_19', 'yjy_17_f_23', 'lym_18_f_19', 'pr_19_f_23', 'ws_20_f_21', ...
                 'wn_21_f_21', 'hjy_22_f_18', 'qyk_23_f_22', 'yd_24_f_20'};
    subjLab   = {'wsn1', 'dy2', 'haq3', 'hry4', 'zjx5', 'yyq6', 'zkw7', 'zy8', 'hys9', 'cjj10', ...
                 'dwq11', 'ljl12', 'jyx13', 'zk14', 'lsy15', 'cjl16', 'yjy17', 'lym18', 'pr19', 'ws20', ...
                 'wn21', 'hjy22', 'qyk23', 'yd24'};
    subjAges    = [18, 22, 24, 20, 20, 18, 18, 20, 20, 18, ...
                   22, 20, 19, 21, 19, 19, 23, 19, 23, 21, ...
                   21, 18, 22, 20];
    subjGenders = {'F', 'F', 'F', 'F', 'M', 'F', 'M', 'F', 'M', 'M', ...
                   'F', 'M', 'M', 'F', 'M', 'M', 'F', 'F', 'F', 'F', ...
                   'F', 'F', 'F', 'F'};

elseif isequal(ExpWord, 'ExplicitExp')
    subj_list = {'hsp_1_m_21', 'hr_2_f_22', 'pxy_3_m_19', 'wyx_4_m_22', 'lml_5_f_19', 'lf_6_m_21', 'md_7_f_21', 'srz_8_m_20', 'fsq_9_f_18', 'caq_10_f_20', ...
                 'xjm_11_f_26', 'cjs_12_m_19', 'xxx_13_f_19', 'wsq_14_f_19', 'zzm_15_f_18', 'lsy_16_f_20', 'man_17_f_19', 'zxy_18_f_18', 'cyh_19_f_21', 'szn_20_f_19', ...
                 'sb_21_m_22', 'drq_22_f_18', 'cy_23_f_18', 'zhc_24_m_21'}; 
    subjLab   = {'hsp1', 'hr2', 'pxy3', 'wyx4', 'lml5', 'lf6', 'md7', 'srz8', 'fsq9', 'caq10', ...
                 'xjm11', 'cjs12', 'xxx13', 'wsq14', 'zzm15', 'lsy16', 'man17', 'zxy18', 'cyh19', 'szn20', ...
                 'sb21', 'drq22', 'cy23', 'zhc24'};  
    subjAges    = [21, 22, 19, 22, 19, 21, 21, 20, 18, 20, ...
                   26, 19, 19, 19, 18, 20, 19, 18, 21, 19, ...
                   22, 18, 18, 21];
    subjGenders = {'M', 'F', 'M', 'M', 'F', 'M', 'F', 'M', 'F', 'F', ...
                   'F', 'M', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', ...
                   'M', 'F', 'F', 'M'};

elseif isequal(ExpWord, 'ImplicitRandExp')
    subj_list = {'zyh_1_f_21', 'why_2_m_18', 'cr_3_f_19', 'zyx_4_f_25', 'wym_5_m_23', 'wd_6_f_18', 'lyh_7_f_21', 'zr_8_f_24', 'zyh_9_m_19', 'zzy_10_f_19', ...
                 'smq_11_f_25', 'sz_12_f_20', 'lzy_13_m_22', 'yxy_14_f_18', 'zxl_15_f_18', 'wqh_16_m_23', 'zxj_17_m_21', 'skx_18_f_20', 'zlh_19_m_24', 'gwt_20_f_23', ...
                 'lwn_21_f_18', 'lrp_22_m_21', 'sjj_23_f_18', 'xy_24_f_19'}; % , 'xr_18_f_23'
    subjLab   = {'zyh1', 'why2', 'cr3', 'zyx4', 'wym5', 'wd6', 'lyh7', 'zr8', 'lyh9', 'zzy10', ...
                 'smq11', 'sz12', 'lzy13', 'yxy14', 'zxl15', 'wqh16', 'zxj17', 'skx18', 'zlh19', 'gwt20', ...
                 'lwn21', 'lrp22', 'sjj23', 'xy24'}; % , 'xr18'
    subjAges    = [21, 18, 19, 25, 23, 18, 21, 24, 19, 19, ...
                   25, 20, 22, 18, 18, 23, 21, 20, 24, 23, ...
                   18, 21, 18, 19];
    subjGenders = {'F', 'M', 'F', 'F', 'M', 'F', 'F', 'F', 'M', 'F', ...
                   'F', 'F', 'M', 'F', 'F', 'M', 'M', 'F', 'M', 'F', ...
                   'F', 'M', 'F', 'F'};
end
subjNum   = 1 : length(subj_list);
subj_listBv = subj_list;
subLen      = length(subj_list);

%% quantify the age and gender
disp('---------- Age of YAs ----------');
[mean_age, sem_age] = Mean_and_Se(subjAges')

disp('---------- Females in YAs----------');
sum(strcmp(subjGenders, 'F'))
length(subjGenders)


%% Experiment parameters
folder = '/Users/ren/Projects-NeuroCode/MyExperiment/HierarchicalCluster';
if isequal(ExpWord, 'ImplicitExp')
    bhvDataDir = [folder, '/FormalExp-LynnNetwork-Results/']; 
elseif isequal(ExpWord, 'ExplicitExp')
    bhvDataDir = [folder, '/FormalExp-ExplicitLearning-LynnNetwork-Results/'];
elseif isequal(ExpWord, 'ImplicitRandExp')
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
predRatio = zeros(subLen, nBlock);
errRatio  = zeros(subLen, nBlock); %% ratio of error or time out trials
fastRatio = zeros(subLen, nBlock);
corrPreRt = zeros(subLen, 1); %% correction ratio of prediction trials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
baseAcc      = zeros(subLen, 1);
transTpAcc   = zeros(subLen, 2); %% correct response ratio in within- and between-cluster transitions
transTpAcc_lure = zeros(subLen, 2, 2); % the 1st and 2nd '2' denote 'lure distractor exists vs. none' and 'transition from boundary node to within node vs. from boundary to boundary'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RTsubj              = zeros(nTrials, subLen);
RTtransTp           = zeros(subLen, 2);
RTtransTp_type      = zeros(subLen, 2, 2); % the 1st and 2nd '2' denote 'within vs. between cluster transition' and 'random and hamiltonian' walk
RTtransTp_Blc       = zeros(subLen, 2, nBlock);
RTtransTp_type_Blc  = cell(1, 2); % random vs. hamiltonian walks
RTtransTp_lure      = zeros(subLen, 2, 2); % the 1st and 2nd '2' denote 'lure distractor exists vs. none' and 'transition from boundary node to within node vs. from boundary to boundary'
lureTrialLen        = zeros(subLen, 2 * 2);
% first 2:  with or without lure stimulus
% second 2: boundary-to-within transition or bounary-to-boundary transition
% third 3:  random or hamiltonian walk
RTtransTp_type_lure = zeros(subLen, 2, 2, 2);
RToverall_Blc_dtNum = zeros(subLen, 3, nBlock); % 3 distractor conditions
lme_subj = [];
for SubIdx = 1 : subLen
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
    corrIdx      = find(errorIdx_Col ~= 1);
    respRT_tmp   = respRT_Col;
    respRT_tmp(errorIdx_Col == 1) = nan;
    RTsubj(:, SubIdx) = respRT_tmp;
    %%% ratio of correct responses through predictions acorss 7 blocks
    for iBlock = 1 : nBlock
        ttBlcIdx  = (trialsInBlc(iBlock, 1) : trialsInBlc(iBlock, 2))';
        trialsBlc = length(ttBlcIdx);
        if isequal(ExpWord, 'ImplicitExp') || isequal(ExpWord, 'ImplicitRandExp')
            predRatio(SubIdx, iBlock) = length(find(respRT_Col(ttBlcIdx) <= grayoffT & errorIdx_Col(ttBlcIdx) ~= 1 & fast_Col(ttBlcIdx) ~= 1)) / trialsBlc; %% equaling to the ratio of 100 bonus
        elseif isequal(ExpWord, 'ExplicitExp')
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
    % 2nd dimension: lure stimulus exists or not;
    % 3rd dimension: within-cluster transition or between-cluster
    % transition
    transTpAcc_lure(SubIdx, 1, 1) = length(find(nodesLabel == 1 & transStyle == 1 & lureIn == 1 & errorIdx_Col ~= 1 & fast_Col ~= 1)) / length(find(nodesLabel == 1 & transStyle == 1 & lureIn == 1));
    transTpAcc_lure(SubIdx, 1, 2) = length(find(nodesLabel == 1 & transStyle ~= 1 & lureIn == 1 & errorIdx_Col ~= 1 & fast_Col ~= 1)) / length(find(nodesLabel == 1 & transStyle ~= 1 & lureIn == 1));
    transTpAcc_lure(SubIdx, 2, 1) = length(find(nodesLabel == 1 & transStyle == 1 & lureIn == 0 & errorIdx_Col ~= 1 & fast_Col ~= 1)) / length(find(nodesLabel == 1 & transStyle == 1 & lureIn == 0));
    transTpAcc_lure(SubIdx, 2, 2) = length(find(nodesLabel == 1 & transStyle ~= 1 & lureIn == 0 & errorIdx_Col ~= 1 & fast_Col ~= 1)) / length(find(nodesLabel == 1 & transStyle ~= 1 & lureIn == 0));
    
    %% RT Comparison: lure distractor exists vs. none
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
    % ------ lure stimulus exists ------
    RTtransTp_type_lure(SubIdx, 1, 1, 1) = mean(respRT_Col(find(nodesLabel == 1 & transStyle == 1 & lureIn == 1 & errorIdx_Col ~= 1 & fast_Col ~= 1 & rndHam_Col == 1))); %% rndHam_Col == 1: random trial
    RTtransTp_type_lure(SubIdx, 1, 2, 1) = mean(respRT_Col(find(nodesLabel == 1 & transStyle ~= 1 & lureIn == 1 & errorIdx_Col ~= 1 & fast_Col ~= 1 & rndHam_Col == 1))); 
    RTtransTp_type_lure(SubIdx, 1, 1, 2) = mean(respRT_Col(find(nodesLabel == 1 & transStyle == 1 & lureIn == 1 & errorIdx_Col ~= 1 & fast_Col ~= 1 & rndHam_Col == 2))); %% rndHam_Col == 2: hamiltonian trial
    RTtransTp_type_lure(SubIdx, 1, 2, 2) = mean(respRT_Col(find(nodesLabel == 1 & transStyle ~= 1 & lureIn == 1 & errorIdx_Col ~= 1 & fast_Col ~= 1 & rndHam_Col == 2))); 
    % ------ no lure stimulus ------
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
    %% within vs. between cluster transition RTs across blocks
    for iBlock = 1 : nBlock
        ttBlcIdx = (trialsInBlc(iBlock, 1) : trialsInBlc(iBlock, 2))';
        RTtransTp_Blc(SubIdx, 1, iBlock) = mean(respRT_Col(find((trials_Col >= ttBlcIdx(1) & trials_Col <= ttBlcIdx(end)) & transStyle == 1 & errorIdx_Col ~= 1 & fast_Col ~= 1)));
        RTtransTp_Blc(SubIdx, 2, iBlock) = mean(respRT_Col(find((trials_Col >= ttBlcIdx(1) & trials_Col <= ttBlcIdx(end)) & transStyle ~= 1 & errorIdx_Col ~= 1 & fast_Col ~= 1)));
    end
    %% within vs. between cluster transition RTs across blocks seperately for random and hamiltonian walk
    for rOh = 1 : 2 % RTtransTp_type_Blc = cell(1, 2);
        for iBlock = 1 : nBlock
            ttBlcIdx = (trialsInBlc(iBlock, 1) : trialsInBlc(iBlock, 2))';
            RTtransTp_type_Blc{rOh}(SubIdx, 1, iBlock) = mean(respRT_Col(find((trials_Col >= ttBlcIdx(1) & trials_Col <= ttBlcIdx(end)) & transStyle == 1 & errorIdx_Col ~= 1 & fast_Col ~= 1 & rndHam_Col == rOh)));
            RTtransTp_type_Blc{rOh}(SubIdx, 2, iBlock) = mean(respRT_Col(find((trials_Col >= ttBlcIdx(1) & trials_Col <= ttBlcIdx(end)) & transStyle ~= 1 & errorIdx_Col ~= 1 & fast_Col ~= 1 & rndHam_Col == rOh)));
        end
    end
    %% RTs with different distractor numbers
    dtNum = dtNum - 1;
    for iDt = 1 : length(dtrNums) % RToverall_Blc_dtNum = zeros(subLen, 3, nBlock); % 3 distractor conditions
        for iBlock = 1 : nBlock
            ttBlcIdx = (trialsInBlc(iBlock, 1) : trialsInBlc(iBlock, 2))';
            RToverall_Blc_dtNum(SubIdx, iDt, iBlock) = mean(respRT_Col(find((trials_Col >= ttBlcIdx(1) & trials_Col <= ttBlcIdx(end)) & dtNum == dtrNums(iDt) & errorIdx_Col ~= 1 & fast_Col ~= 1)));
        end
    end

    %% table for the lme test: RT, 
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
    lme_i = [repmat(SubIdx, nTrials, 1), trials_Col, respRT_Col, from_nodes, to_nodes, dtNum, is_BdNode, recency_nodes, transStyle, rndHam_Col, lureIn];
    lme_i = lme_i(corrIdx, :);
    lme_i(:, 2) = log(lme_i(:, 2)); %% log(trial)
    lme_subj = [lme_subj; lme_i];
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


%% SI figure, RT Learning effect: within- and between-trans in Random and Hamiltonian Walk
% RTtransTp_type_Blc  = cell(1, 2); % random vs. hamiltonian walks
% RTtransTp_type_Blc{1} = (subLen, transType, nBlock)
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
figure('Position', [100 100 260 120]), clf;
hold on;
nWalks = 2; % Random and Hamiltonian Walk
nTrans = 2; % within and Between-cluster transition
for i_rOh = 1 : nWalks
    RTtrans_iWs = RTtransTp_type_Blc{i_rOh};
    for iTs = 1 : nTrans
        RTtrans_iWs_iTs = squeeze(RTtrans_iWs(:, iTs, :));
        [RTsAvg, RTsSem] = Mean_and_Se(RTtrans_iWs_iTs, 1);
        if i_rOh == 1
            colorTmp = colorSet(iTs, :);
        elseif i_rOh == 2
            colorTmp = 0.5 * colorSet(iTs, :) + 0.5 * [1, 1, 1];
        end
        errorbar(1 : 1 : nBlock, RTsAvg, RTsSem, 'Color', colorTmp, 'LineStyle', '-', 'LineWidth', errLineWid); hold on;
        for iBlc = 1 : nBlock
            plot(iBlc, RTsAvg(iBlc), 'Marker', 'o', 'MarkerSize', markSize, 'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', colorTmp, 'LineStyle', '-'); hold on;
        end
    end
end
ylim([1, 1.6]);
xlim([0.5, nBlock+0.5]);
if figKey == 0
    % ------For presentation------
    set(gca, 'LineWidth', 2);
    set(gca, 'FontSize', 15, 'FontWeight', 'bold', 'FontName', 'Arial');
    set(gca, 'XTick', 1 : 1 : nBlock, 'XTickLabel', 1 : 1 : nBlock);
    set(gca, 'YTick', [1, 1.3, 1.6], 'YTickLabel', [1, 1.3, 1.6]); 
elseif figKey == 1
    % ------For Adobe Illustrator------
    set(gca, 'LineWidth', 0.8);
    set(gca, 'FontSize', 10, 'FontWeight', 'bold', 'FontName', 'Arial');
    set(gca, 'XTick', 1 : 1 : nBlock, 'XTickLabel', '');
    set(gca, 'YTick', [1, 1.3, 1.6], 'YTickLabel', {'', '', ''}); 
end
box off;

%% SI figure, RT Learning effect: different distractor numbers
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
figure('Position', [100 100 260 120]), clf;
hold on;
LineStys = {'-', '-.', ':'};
for iDt = 1 : length(dtrNums) % RToverall_Blc_dtNum = zeros(subLen, 3, nBlock); % 3 distractor conditions
    RTtrans_iDt = squeeze(RToverall_Blc_dtNum(:, iDt, :));
    [RTsAvg, RTsSem] = Mean_and_Se(RTtrans_iDt, 1);
    errorbar(1 : 1 : nBlock, RTsAvg, RTsSem, 'Color', blueGrad(iDt, :), 'LineStyle', LineStys{iDt}, 'LineWidth', errLineWid); hold on;
    for iBlc = 1 : nBlock
        plot(iBlc, RTsAvg(iBlc), 'Marker', 'o', 'MarkerSize', markSize, 'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', blueGrad(iDt, :), 'LineStyle', '-'); hold on;
    end
end
xlim([0.5, nBlock+0.5]);
ylim([0.9, 1.6]);
if figKey == 0
    % ------For presentation------
    set(gca, 'LineWidth', 2);
    set(gca, 'FontSize', 15, 'FontWeight', 'bold', 'FontName', 'Arial');
    set(gca, 'XTick', 1 : 1 : nBlock, 'XTickLabel', 1 : 1 : nBlock);
    set(gca, 'YTick', [1, 1.3, 1.6], 'YTickLabel', [1, 1.3, 1.6]); 
elseif figKey == 1
    % ------For Adobe Illustrator------
    set(gca, 'LineWidth', 0.8);
    set(gca, 'FontSize', 10, 'FontWeight', 'bold', 'FontName', 'Arial');
    set(gca, 'XTick', 1 : 1 : nBlock, 'XTickLabel', '');
    set(gca, 'YTick', [0.9, 1, 1.3, 1.6], 'YTickLabel', {'', '', '', ''}); 
end
box off;

%% BehavioralPaper, Figure xx: overall RTs for within vs. between transitions in Random and Hamiltonian Walk
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
figure('Position', [100 100 260 120]), clf;

% RTtransTp_type = zeros(subLen, 2, 2); % the 1st and 2nd '2' denote 'within vs. between cluster transition' and 'random and hamiltonian' walk
[RTavg_rOh, RTsem_rOh] = Mean_and_Se(RTtransTp_type, 1);
RTavg_rOh = squeeze(RTavg_rOh);
RTsem_rOh = squeeze(RTsem_rOh);
barPos = [1, 1.5; 1.7, 2.2];
for i_rOh = 1 : 2 % Random vs. Hamiltonian Walk
    if i_rOh == 1
        walkWord = 'Radom';
    elseif i_rOh == 2
        walkWord = 'Hamiltonian';
    end
    barPos_i = barPos(i_rOh, :);
    %%% line plot
    plot(barPos_i, RTtransTp_type(:, :, i_rOh), 'Color', [0.6, 0.6, 0.6], 'LineStyle', '-', 'LineWidth', indvLineW); hold on;
    plot(barPos_i, RTavg_rOh(:, i_rOh), 'Color', [0, 0, 0], 'LineStyle', '-', 'LineWidth', errLineWid); hold on;
    for i_tp = 1 : 2 % within vs. between-transition
        if i_rOh == 1
            colorTmp = colorSet(i_tp, :);
        elseif i_rOh == 2
            colorTmp = 0.5 * colorSet(i_tp, :) + 0.5 * [1, 1, 1];
        end
        errorbar(barPos_i(i_tp), RTavg_rOh(i_tp, i_rOh), RTsem_rOh(i_tp, i_rOh), 'Color', 'k', 'LineStyle', 'none', 'LineWidth', errLineWid); hold on;
        plot(barPos_i(i_tp), RTavg_rOh(i_tp, i_rOh), 'Marker', 'o', 'MarkerSize', markSize, 'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', colorTmp, 'LineStyle', '-'); hold on;
    end
    % ------ Statistical tests ------
    disp(['======== ', walkWord, '-within vs. between trans ========']);
    [h, p, ci, stats] = ttest(RTtransTp_type(:, 1, i_rOh), RTtransTp_type(:, 2, i_rOh))
end
xlim([0.6, 2.6]);
ylim([0.5, 2]);
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
    set(gca, 'YTick', [0.5, 1, 1.5, 2], 'YTickLabel', '');
end
box off;

%% the difference: report in the manuscript
[RTavg_rOhDif, RTsem_rOhDif] = Mean_and_Se(RTtransTp_type(:, 2, :) - RTtransTp_type(:, 1, :), 1);
RTavg_rOhDif = squeeze(RTavg_rOhDif);
RTsem_rOhDif = squeeze(RTsem_rOhDif);
figure('Position', [100 100 180 120]), clf;
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
% ---------- Comparison between within vs. between-trans within a walk type ----------
disp('---------- Random Walk: within vs. between-trans ----------');
[h_r, p_r, ci_r, stats_r] = ttest(RTtransTp_type(:, 2, 1), RTtransTp_type(:, 1, 1))
disp('---------- Hamiltonian Walk: within vs. between-trans ----------');
[h_h, p_h, ci_h, stats_h] = ttest(RTtransTp_type(:, 2, 2), RTtransTp_type(:, 1, 2))
% ---------- whether there is significant difference between random and hamiltonian walk ----------
disp('---------- Random vs. Hamiltonian Walk: Diff in between vs. within trans ----------');
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

%% with lure distractors vs. without
% RTtransTp_lure      = zeros(subLen, 2, 2); % the 1st and 2nd '2' denote 'lure distractor exists vs. none' and 'transition from boundary node to within node vs. from boundary to boundary'
% RTtransTp_type_lure = zeros(subLen, 2, 2, 2); % the last 2 denotes 'random and hamiltonian' walk
flg = 3;
if flg == 1
    RTmat = RTtransTp_lure; %% mixture of random and hamiltonian path
elseif flg == 2
    RTmat = RTtransTp_type_lure(:, :, :, 1); %% only random path
elseif flg == 3
    RTmat = RTtransTp_type_lure(:, :, :, 2); %% only hamiltonian path
end
%% BehavioralPaper, Figure xx: influence of lure stimulus on RTs
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
figure('Position', [100 100 260 120]), clf;

[RTavg_rOh, RTsem_rOh] = Mean_and_Se(RTmat, 1);
RTavg_rOh = squeeze(RTavg_rOh);
RTsem_rOh = squeeze(RTsem_rOh);
barPos = [1, 1.5; 1.7, 2.2];
for iTb = 1 : 2     % 'transition from boundary node to within node vs. from boundary to boundary'
    if iTb == 1
        transWord = 'boundary-to-within';
    elseif iTb == 2
        transWord = 'boundary-to-boundary';
    end
    barPos_i = barPos(iTb, :);
    %%% line plot
    plot(barPos_i, RTmat(:, :, iTb), 'Color', [0.6, 0.6, 0.6], 'LineStyle', '-', 'LineWidth', indvLineW); hold on;
    plot(barPos_i, RTavg_rOh(:, iTb), 'Color', [0, 0, 0], 'LineStyle', '-', 'LineWidth', errLineWid); hold on;
    for ilr = 1 : 2 % 'lure distractor exists vs. none'
        if ilr == 1
            colorTmp = [0, 0, 0]; % with lure distractor
        elseif ilr == 2
            colorTmp = [1, 1, 1]; % without lure distractor
        end
        errorbar(barPos_i(ilr), RTavg_rOh(ilr, iTb), RTsem_rOh(ilr, iTb), 'Color', 'k', 'LineStyle', 'none', 'LineWidth', errLineWid); hold on;
        plot(barPos_i(ilr), RTavg_rOh(ilr, iTb), 'Marker', 'o', 'MarkerSize', markSize, 'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', colorTmp, 'LineStyle', '-'); hold on;
    end
    % ------ Statistical tests ------
    disp(['======== ', transWord, ': with vs. without lure stimulus ========']);
    [h, p, ci, stats] = ttest(RTmat(:, 1, iTb), RTmat(:, 2, iTb))
end
xlim([0.6, 2.6]);
ylim([0.5, 2]);
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
    set(gca, 'YTick', [0.5, 1, 1.5, 2], 'YTickLabel', '');
end
box off;


%% Quantify the RT difference when the lure stimulus exists vs. not
[RTavg_rOhDif, RTsem_rOhDif] = Mean_and_Se(RTmat(:, 2, :) - RTmat(:, 1, :), 1);
RTavg_rOhDif = squeeze(RTavg_rOhDif);
RTsem_rOhDif = squeeze(RTsem_rOhDif);
figure('Position', [100 100 180 120]), clf;
for iB = 1 : 2
    if iB == 1     %% transition from boundary node to within node
        LineBar = '-';
    elseif iB == 2 %% from boundary to boundary
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
% ---------- Comparison between within vs. between-trans within a walk type ----------
disp('---------- Boundary-to-within: with vs. without lure ----------');
[h_r, p_r, ci_r, stats_r] = ttest(RTmat(:, 2, 1), RTmat(:, 1, 1))
disp('---------- Boundary-to-boundary: with vs. without lure ----------');
[h_h, p_h, ci_h, stats_h] = ttest(RTmat(:, 2, 2), RTmat(:, 1, 2))
% ---------- whether there is significant difference between random and hamiltonian walk ----------
disp('---------- Boundary-to-within vs. boundar-to-boundary: Diff in with vs. without lure ----------');
[h_rOh, p_rOh, ci_rOh, stats_rOh] = ttest(RTmat(:, 2, 2)-RTmat(:, 1, 2), RTmat(:, 2, 1)-RTmat(:, 1, 1))


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
%     save([subjDir, subjBv, 'mouseTraj_trials_blc.mat'], 'mouseTraj_trials_blc');
%     save([subjDir, subjBv, 'respYes_trials_blc.mat'], 'respYes_trials_blc');
end









