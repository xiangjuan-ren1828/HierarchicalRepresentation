% ClusterRep_v3_lureDistractor
% write by rxj @ 01/09/2020 
% In some trials, a lure distractor will be presented simultaneously with
% the target, both of which belong to the same cluster but without
% connections

%%
clear
clc

%%
cd ../
addpath('tight_subplot/');
cd 'HierarchicalCluster/'

addpath('clusterRep-v2-dataAnal/');
addpath('clusterRep-v3-dataAnal/');
addpath('clusterRep-version4-spatial/');

%% Subject
ExpWord_List = {'cueEyeTrack1', 'cueEyeTrack2', 'cueEyeTrack3', 'cueEyeTrack4', 'cueEyeTrack5'};
ExpIdx       = 5;
ExpWord      = ExpWord_List{ExpIdx};
if isequal(ExpWord, 'cueEyeTrack1')
    subj_list = {'yjw_2_f_21', 'lyc_3_m_20', 'lmh_5_f_19', 'lh_6_f_19', 'stj_7_f_20', 'ljl_8_m_21'}; % No eyetracking data: 'mcaj_1_f_20', 'cjq_4_f_18', 
    subjLab   = {'yjw2', 'lyc3', 'lmh5', 'lh6', 'stj7', 'ljl8'}; % 'mcaj1', 'cjq4', 
elseif isequal(ExpWord, 'cueEyeTrack2')
    subj_list = {'csx_10_m_20', 'ljj_11_f_24', 'lx_12_m_19', 'jzc_14_m_24', 'mym_17_f_18'};  
    subjLab   = {'csx10', 'ljj11', 'lx12', 'jzc14', 'mym17'};
elseif isequal(ExpWord, 'cueEyeTrack3')
    subj_list = {'zxh_22_f_20', 'ljy_23_f_21', 'sq_24_f_25'};  
    subjLab   = {'zxh22', 'ljy23', 'sq24'};
elseif isequal(ExpWord, 'cueEyeTrack4')
    subj_list = {'ftt_25_f_21', 'wyg_26_f_20', 'lx_27_f_21', 'wgj_28_m_19', 'jj_29_m_18', 'rxj_1_f_01'};
    subjLab   = {'ftt25', 'wyg26', 'lx27', 'wgj28', 'jj29' , 'rxj1'}; 
elseif isequal(ExpWord, 'cueEyeTrack5')
    subj_list = {'lxy_1_f_23', 'hyc_2_m_18', 'zzx_3_m_23', 'zjy_4_f_20', 'lyw_5_f_24', 'cx_6_m_22'}; % 'rxj_2_f_01',  
    subjLab   = {'lxy1', 'hyc2', 'zzx3', 'zjy4', 'lyw5', 'cx6'}; % 'rxj2', 
end
subjNum   = 1 : length(subj_list);
subj_listBv = subj_list;
subLen      = length(subj_list);

%% Experiment parameters
folder     = pwd;
bhvDataDir = [folder, '/EyeTrackFormalExp-LynnNetwork-Results/'];
eyeDataDir = [folder, '/EyeTrackFormalExp-LynnNetwork-EyeTrackResults/'];
rndTrial = 700;
HamTrial = 800; 
nTrials  = rndTrial + HamTrial;
rndInter = 85;
HamInter = 15;

tInBlcMouse = [230, 2*230, 2*230+240, rndTrial + (200 : 200 : HamTrial)];
tInBlcEye   = [nTrials + (200 : 200 : HamTrial)];
nBlockMouse = length(tInBlcMouse); %% seperate the total trials into 10 blocks
nBlockEye   = length(tInBlcEye);
nBlock      = nBlockMouse + nBlockEye;
tInBlc      = [tInBlcMouse, tInBlcEye];
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

% Durations
if isequal(ExpWord, 'cueEyeTrack1')
    fixationDur = 0.6; %% seconds
    respDur     = 2;   %3s;
    grayoffT    = 0.8; %1s;   %% For the 1st second, the gray cue won't appear.
    grayOnT     = 0.8; %1s;   %% For 1-2s, the gray square will be overlapped on the target image; for 2-3s, the gray square superimposed on the target image with the largest gray scale
    bonusDec    = 0.5; %% bonus decreased from 100 to 0 within 0.5s
    bonusDecEye = 0.2;
    bonusTime   = 1;   %% Only for practice trials
    saccTole    = 0.4;
elseif isequal(ExpWord, 'cueEyeTrack2') || isequal(ExpWord, 'cueEyeTrack3') || isequal(ExpWord, 'cueEyeTrack4')
    fixationDur = 0.6; %% seconds
    respDur     = 2;   %3s;
    grayoffT    = 1; %1s;   %% For the 1st second, the gray cue won't appear.
    grayOnT     = 0.8; %1s;   %% For 1-2s, the gray square will be overlapped on the target image; for 2-3s, the gray square superimposed on the target image with the largest gray scale
    bonusDec    = 0.5; %% bonus decreased from 100 to 0 within 0.5s
    bonusDecEye = 0.2;
    bonusTime   = 1;   %% Only for practice trials
    saccTole    = 0.3;
elseif isequal(ExpWord, 'cueEyeTrack5')
    fixationDur = 0.6; %% seconds
    respDur     = 2.5;   %3s;
    grayoffT    = 1; %1s;   %% For the 1st second, the gray cue won't appear.
    grayOnT     = 0.8; %1s;   %% For 1-2s, the gray square will be overlapped on the target image; for 2-3s, the gray square superimposed on the target image with the largest gray scale
    bonusDec    = 0.5; %% bonus decreased from 100 to 0 within 0.5s
    bonusDecEye = 0.3;
    bonusTime   = 1;   %% Only for practice trials
    saccTole    = 0.5;
end

%% Columns index for the clusterResult.mat
if isequal(ExpWord, 'cueEyeTrack1') || isequal(ExpWord, 'cueEyeTrack2') || isequal(ExpWord, 'cueEyeTrack3') || isequal(ExpWord, 'cueEyeTrack4')
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
    imgAngs  = 26 : 29;
    fixLast  = 30;
    base_err = 31 : 32;
    resp_err = 33 : 34;
    fixFroze = 35;
elseif isequal(ExpWord, 'cueEyeTrack5')
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
    errCnt   = 15;
    bonusCn  = 16; %% bonus in each trial
    trialLt  = 17; %% trial total length
    imgPosX  = 18 : 21;
    imgPosY  = 22 : 25;
    fastIdx  = 26;
    imgAngs  = 27 : 30;
    fixLast  = 31;
    base_err = 32 : 33;
    resp_err = 34 : 35;
    fixFroze = 36;
end

%% Analyzing the data
if isequal(ExpWord, 'cueEyeTrack1') || isequal(ExpWord, 'cueEyeTrack2') || isequal(ExpWord, 'cueEyeTrack3')
    TotBlocks = nBlockEye;
    TotTrials = HamTrial;
elseif isequal(ExpWord, 'cueEyeTrack4') || isequal(ExpWord, 'cueEyeTrack5')
    TotBlocks = nBlock;
    TotTrials = nTrials + HamTrial;
end
RTtransTp_lure = zeros(subLen, 2, 2); % the 1st and 2nd '2' denote 'lure distractor exists vs. none' and 'transition from boundary node to within node vs. from boundary to boundary'
lureTrialLen   = zeros(subLen, 2 * 2);
RTtransTp_type_lure = zeros(subLen, 2, 2, 2); % the last 2 denotes 'random and hamiltonian' walk
lme_subj = [];

trlInB = 100; %% trial No. in each bins
nBins  = TotTrials / trlInB;
RTtransTp_bin      = nan(subLen, 2, nBins);
RTtransTp_type_bin = cell(1, 2); % random vs. hamiltonian walks
saccRT = saccTole;
for SubIdx = 1 : subLen
    subjBv      = subj_listBv{SubIdx};
    subTit      = subjLab{SubIdx};
    subjDir     = [bhvDataDir, subjBv, '/'];
    if isequal(subTit, 'hyc2') || isequal(subTit, 'zzx3')
        if isequal(subTit, 'hyc2')
            TotBlocks = 9;
        elseif isequal(subTit, 'zzx3')
            TotBlocks = 7;
        end
        clusterResult = [];
        for iCb = 1 : TotBlocks
            load([subjDir, subjBv, 'clusterResult_trials_in_Blc', num2str(iCb), '.mat'], 'clusterResult_trials');
            clusterResult = [clusterResult; clusterResult_trials];
        end
        trialCate     = (1 : 1 : size(clusterResult, 1))';
        clusterResult = [trialCate, clusterResult];
    else
        TotBlocks = nBlock;
        load([subjDir, subjBv, 'clusterResult_BlocksEyeTrack.mat'], 'clusterResult');
    end
    TotTrials = size(clusterResult, 1);
    nBins     = TotTrials / trlInB;
    
    trials_Col   = (1 : 1 : size(clusterResult, 1))';
    respRT_Col   = clusterResult(trials_Col, respRT);
    errorIdx_Col = clusterResult(trials_Col, errorIdx);
    timeOIdx_Col = clusterResult(trials_Col, timeOIdx);
    bonus_Col    = clusterResult(trials_Col, bonusCn);
    rndHam_Col   = clusterResult(trials_Col, trialTp); 
    fast_Col     = clusterResult(trials_Col, fastIdx);
    fixFail_Col  = clusterResult(trials_Col, fixFroze);
    rndHam_Col(rndHam_Col == 1) = 0; %% random trial
    rndHam_Col(rndHam_Col == 2) = 1; %% hamiltonian trial
    corrIdx      = find(errorIdx_Col ~= 1);
  
    %%% within vs. between cluster transition RTs
    from_nodes = clusterResult(trials_Col, objTgt);
    to_nodes   = clusterResult(trials_Col, objResp(1));
    dt_nodes   = clusterResult(trials_Col, objResp);
    transIn    = zeros(TotTrials, 3);
    lureIn     = zeros(TotTrials, 1); % whether there are two items selected from the same cluster
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
    nodesLabel = zeros(TotTrials, 1);
    nodesLabel(from_nodes == 1 | from_nodes == 5 | from_nodes == 6 | from_nodes == 10 | from_nodes == 11 | from_nodes == 15) = 1;
    dtNum = arrayfun(@(x) length(find(dt_nodes(x, :)~=0)), (1 : TotTrials)');
    is_BdNode = zeros(TotTrials, 1); %% whether the "from" nodes are boundary nodes
    is_BdNode(find(from_nodes == 1 | from_nodes == 5 | from_nodes == 6 | from_nodes == 10 | from_nodes == 11 | from_nodes == 15)) = 1;
    to_BdNode = zeros(TotTrials, 1); %% whether the "to" nodes are boundary nodes
    to_BdNode(find(to_nodes == 1 | to_nodes == 5 | to_nodes == 6 | to_nodes == 10 | to_nodes == 11 | to_nodes == 15)) = 1;
    
    %% Comparison: lure distractor exists vs. none
    dtNos = 4;
    %RTtransTp_lure = zeros(subLen, 2, 2); % the 1st and 2nd '2' denote 'lure distractor exists vs. none' and 'transition from boundary node to within node vs. from boundary to boundary'
    RTtransTp_lure(SubIdx, 1, 1) = mean(respRT_Col(find(nodesLabel == 1 & transStyle == 1 & lureIn == 1 & errorIdx_Col ~= 1))); %  & dtNum == dtNos
    RTtransTp_lure(SubIdx, 1, 2) = mean(respRT_Col(find(nodesLabel == 1 & transStyle ~= 1 & lureIn == 1 & errorIdx_Col ~= 1)));
    RTtransTp_lure(SubIdx, 2, 1) = mean(respRT_Col(find(nodesLabel == 1 & transStyle == 1 & lureIn == 0 & errorIdx_Col ~= 1)));
    RTtransTp_lure(SubIdx, 2, 2) = mean(respRT_Col(find(nodesLabel == 1 & transStyle ~= 1 & lureIn == 0 & errorIdx_Col ~= 1)));
    
    lureTrialLen(SubIdx, 1) = length(find(nodesLabel == 1 & transStyle == 1 & lureIn == 1 & errorIdx_Col ~= 1));
    lureTrialLen(SubIdx, 2) = length(find(nodesLabel == 1 & transStyle ~= 1 & lureIn == 1 & errorIdx_Col ~= 1));
    lureTrialLen(SubIdx, 3) = length(find(nodesLabel == 1 & transStyle == 1 & lureIn == 0 & errorIdx_Col ~= 1));
    lureTrialLen(SubIdx, 4) = length(find(nodesLabel == 1 & transStyle ~= 1 & lureIn == 0 & errorIdx_Col ~= 1));
    
    %% Comparison: lure distractor exists vs. none, in random and hamiltonian path
    %RTtransTp_type_lure = zeros(subLen, 2, 2, 2); % the last 2 denotes 'random and hamiltonian' walk
    %%% wtihin vs. between cluster transition RTs separately for random and
    %%% hamiltonian walk
    RTtransTp_type_lure(SubIdx, 1, 1, 1) = mean(respRT_Col(find(nodesLabel == 1 & transStyle == 1 & lureIn == 1 & errorIdx_Col ~= 1 & rndHam_Col == 0))); %% rndHam_Col == 1: random trial
    RTtransTp_type_lure(SubIdx, 1, 2, 1) = mean(respRT_Col(find(nodesLabel == 1 & transStyle ~= 1 & lureIn == 1 & errorIdx_Col ~= 1 & rndHam_Col == 0))); 
    RTtransTp_type_lure(SubIdx, 1, 1, 2) = mean(respRT_Col(find(nodesLabel == 1 & transStyle == 1 & lureIn == 1 & errorIdx_Col ~= 1 & rndHam_Col == 1))); %% rndHam_Col == 2: hamiltonian trial
    RTtransTp_type_lure(SubIdx, 1, 2, 2) = mean(respRT_Col(find(nodesLabel == 1 & transStyle ~= 1 & lureIn == 1 & errorIdx_Col ~= 1 & rndHam_Col == 1))); 
    
    RTtransTp_type_lure(SubIdx, 2, 1, 1) = mean(respRT_Col(find(nodesLabel == 1 & transStyle == 1 & lureIn == 0 & errorIdx_Col ~= 1 & rndHam_Col == 0)));
    RTtransTp_type_lure(SubIdx, 2, 2, 1) = mean(respRT_Col(find(nodesLabel == 1 & transStyle ~= 1 & lureIn == 0 & errorIdx_Col ~= 1 & rndHam_Col == 0)));
    RTtransTp_type_lure(SubIdx, 2, 1, 2) = mean(respRT_Col(find(nodesLabel == 1 & transStyle == 1 & lureIn == 0 & errorIdx_Col ~= 1 & rndHam_Col == 1)));
    RTtransTp_type_lure(SubIdx, 2, 2, 2) = mean(respRT_Col(find(nodesLabel == 1 & transStyle ~= 1 & lureIn == 0 & errorIdx_Col ~= 1 & rndHam_Col == 1))); 
    
    %% Linear mixed model
    recency_nodes = zeros(TotTrials, 3);
    for iT = 1 : TotTrials
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
    %%% 4 kinds of codings
    lme_h = [];
    lme_h = [repmat(SubIdx, TotTrials, 1), trials_Col, respRT_Col, dtNum, is_BdNode, to_BdNode, recency_nodes, rndHam_Col, lureIn];
    lme_h = lme_h(corrIdx, :);
    lme_h(:, 2) = log(lme_h(:, 2)); %% log(trial)
    
    lme_subj = [lme_subj; lme_h];
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
nBinsBar = 20;

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

%% statistical testing: LME
lme_lure = lme_subj;
lme_lure = lme_lure(find(lme_lure(:, 5) == 1), :);
lme_lure(:, 5) = [];

colName = {'subID', 'logtrialNo', 'rt', 'dtNum', 'toBdNode', 'recencyL', 'recencySL', 'recencyTL', 'rndORham', 'lureWith'};
RT_tbl = array2table(lme_lure, 'VariableNames', colName);
RT_tbl.toBdNode   = nominal(RT_tbl.toBdNode);
RT_tbl.rndORham   = nominal(RT_tbl.rndORham);
RT_tbl.lureWith   = nominal(RT_tbl.lureWith);

lmeOutput1 = fitlme(RT_tbl, 'rt~logtrialNo + dtNum + recencyL + recencySL + rndORham + lureWith*dtNum + lureWith*toBdNode + (1 + logtrialNo + recencyL + recencySL + rndORham + lureWith*dtNum + lureWith*toBdNode|subID)') %  + logtrial*stage + recency + clusteridx
anova(lmeOutput1, 'DFMethod', 'satterthwaite') % 'DFMethod': 'residual' (default), 'satterthwaite', 'none'A
