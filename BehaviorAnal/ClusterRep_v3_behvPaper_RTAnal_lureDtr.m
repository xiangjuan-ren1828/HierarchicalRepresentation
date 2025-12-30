% ClusterRep_v3_RTAnal_lureDtr
% write by rxj @ 10/25/2021
%
% RT with lure distractors

%%
clear
clc

%%
addpath('tight_subplot/');
addpath(genpath('HierarchicalCluster/'));

%% Subject
ExpWord_List = {'cue', 'explicit', 'cueRand', 'cueEyeTrack'};
ExpIdx       = 3;
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

%% For the boundary node, if another boundary within the same cluster exists, how about the RTs?
% added by rxj @ 08/24/2021
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

%% calculations
dtrCond = 1; % only trials with 1, 2, 3 distractors or 4 (merge all kinds of distractors)
RT_transIn_subj  = zeros(bndNode_Num, 2, subLen);
RT_transOut_subj = zeros(bndNode_Num, 2, subLen);
transIn_Len  = zeros(bndNode_Num, 2, subLen);
transOut_Len = zeros(bndNode_Num, 2, subLen);
for SubIdx = 1 : subLen
    subjBv      = subj_listBv{SubIdx};
    subTit      = subjLab{SubIdx};
    subjDir     = [bhvDataDir, subjBv, '/'];
    load([subjDir, subjBv, 'clusterResult_Blocks.mat'], 'clusterResult');
    
    respRT_Col   = clusterResult(:, respRT);
    errorIdx_Col = clusterResult(:, errorIdx);
    rndHam_Col   = clusterResult(:, trialTp); 
    fast_Col     = clusterResult(:, fastIdx); 
    %%% within vs. between cluster transition RTs
    from_nodes = clusterResult(:, objTgt);
    to_nodes   = clusterResult(:, objResp(1));
    dt_nodes   = clusterResult(:, objResp);
    exId = find(errorIdx_Col ~= 0 | fast_Col == 1);
    %%% within vs. between cluster transition RTs across different distractor numbers
    dtNum = arrayfun(@(x) length(find(dt_nodes(x, :)~=0)), (1 : nTrials)');
    dtNum = dtNum - 1;
    
    respRT_Col(exId, :) = [];
    from_nodes(exId, :) = [];
    to_nodes(exId, :)   = [];
    dt_nodes(exId, :)   = [];
    dtNum(exId, :)      = [];
    if dtrCond == 1 || dtrCond == 2 || dtrCond == 3
        respRT_Col = respRT_Col(dtNum == dtrCond, :);
        from_nodes = from_nodes(dtNum == dtrCond, :);
        to_nodes   = to_nodes(dtNum == dtrCond, :);
        dt_nodes   = dt_nodes(dtNum == dtrCond, 2 : (dtrCond + 1));
    end

    RT_transIn  = zeros(bndNode_Num, 2);
    RT_transOut = zeros(bndNode_Num, 2);
    for i = 1 : bndNode_Num
        %% ---------transition: from boundary node to within node---------
        transIn_InDtr  = transIn_clsInDtr{i}; % transIn_clsInDtr{1} = [1, 2, 5; 1, 3, 5; 1, 4, 5];
        transIn_OutDtr = transIn_clsOutDtr{i}; 
        RTs_InDtr  = [];
        for j = 1 : size(transIn_InDtr, 1)
            pair_j = transIn_InDtr(j, :);
            dt_Yes = sum((dt_nodes == pair_j(3)), 2);
            RTs_InDtr  = [RTs_InDtr; respRT_Col(find(from_nodes == pair_j(1) & to_nodes == pair_j(2) & dt_Yes == 1))];
        end
        RTs_OutDtr = [];
        for k = 1 : size(transIn_OutDtr, 1)
            pair_k = transIn_OutDtr(k, :);
            dt_Yes = sum((dt_nodes == pair_k(3)), 2);
            RTs_OutDtr = [RTs_OutDtr; respRT_Col(find(from_nodes == pair_k(1) & to_nodes == pair_k(2) & dt_Yes == 1))];
        end
        RT_transIn(i, 1) = nanmean(RTs_InDtr);
        RT_transIn(i, 2) = nanmean(RTs_OutDtr);

        %%% trial length of each combination
        transIn_Len(i, 1, SubIdx)  = length(RTs_InDtr);
        transIn_Len(i, 2, SubIdx)  = length(RTs_OutDtr);
    
        %% ---------transition: from boundary node to boundary node---------
        transOut_InDtr  = transOut_clsInDtr{i};
        transOut_OutDtr = transOut_clsOutDtr{i};
        RTs_InDtr = [];
        for j = 1 : size(transOut_InDtr, 1)
            pair_j = transOut_InDtr(j, :);
            dt_Yes = sum((dt_nodes == pair_j(3)), 2);
            RTs_InDtr = [RTs_InDtr; respRT_Col(find(from_nodes == pair_j(1) & to_nodes == pair_j(2) & dt_Yes == 1))];
        end
        
        RTs_OutDtr = [];
        for k = 1 : size(transOut_OutDtr, 1)
            pair_k = transOut_OutDtr(k, :);
            dt_Yes = sum((dt_nodes == pair_k(3)), 2);
            RTs_OutDtr = [RTs_OutDtr; respRT_Col(find(from_nodes == pair_k(1) & to_nodes == pair_k(2) & dt_Yes == 1))];
        end
        RT_transOut(i, 1) = nanmean(RTs_InDtr);
        RT_transOut(i, 2) = nanmean(RTs_OutDtr);
        
        %%% trial length of each combination
        transOut_Len(i, 1, SubIdx)  = length(RTs_InDtr);
        transOut_Len(i, 2, SubIdx)  = length(RTs_OutDtr);
        
    end
    RT_transIn_subj(:, :, SubIdx)  = RT_transIn;
    RT_transOut_subj(:, :, SubIdx) = RT_transOut;
end

%% color settings
colorSet = [249, 183, 176; ... % pink
    248, 218, 172; ...         % orange
    138, 170, 51; ...          % dark green
    210, 234, 200; ...         % light green
    84, 185, 211; ...          % dark blue
    184, 204, 225; ...         % light blue
    198, 127, 192; ...         % dark purple
    219, 204, 226; ...
    102, 102, 102; ...
    255, 255, 255] ./ [255, 255, 255];
colorBehv = colorSet([1, 3], :);
colorSubj = colorSet([1, 2, 3, 9], :);
colorTrans = colorSet([1, 2, 3, 9, 10], :);

% colorFace  = [167, 214, 221; ...
%              61, 142, 143; ...
%              131, 185, 199] ./ [255, 255, 255]; 
% redColor   = 0.6 * [1, 0.64, 0] + 0.4 * [1, 1, 1];


%% average
%%% paired-sample t test
RT_transIn_nodeAvg  = squeeze(nanmean(RT_transIn_subj, 1)); % 2 * subLen
RT_transOut_nodeAvg = squeeze(nanmean(RT_transOut_subj, 1));
RT_transInDif_nodeAvg  = squeeze(nanmean(RT_transIn_subj(:, 2, :) - RT_transIn_subj(:, 1, :), 1)); % 1 * subLen
RT_transOutDif_nodeAvg = squeeze(nanmean(RT_transOut_subj(:, 2, :) - RT_transOut_subj(:, 1, :), 1)); 

%%% group mean
[transIn_avg, transIn_sem]   = Mean_and_Se(squeeze(nanmean(RT_transIn_subj, 1)), 2);
[transOut_avg, transOut_sem] = Mean_and_Se(squeeze(nanmean(RT_transOut_subj, 1)), 2);

%% plotting: within-trans & between-trans separately
% plot in the same figure, and as lines
tail_id = 'both';
barPos = [0.5, 1; 1.2, 1.7];
figure('Position', [100 100 350 200]), clf;
hold on;
for iTrans = 1 : 2 % within- & between-transitions
    barPos_i = barPos(iTrans, :);
    if iTrans == 1
        transAvg = transIn_avg;
        transSem = transIn_sem;
        data_i   = RT_transIn_nodeAvg';
    elseif iTrans == 2
        transAvg = transOut_avg;
        transSem = transOut_sem;
        data_i   = RT_transOut_nodeAvg';
    end
    plot(barPos_i, data_i, '-', 'LineWidth', 2, 'Color', [0.6, 0.6, 0.6]); hold on;
    plot(barPos_i, transAvg, 'Color', colorBehv(iTrans, :), 'Marker', '.', 'LineStyle', '-', 'LineWidth', 4, 'MarkerSize', 1); hold on;
    plot(barPos_i(1), transAvg(1), 'Marker', '.', 'MarkerSize', 10, 'Color', colorBehv(iTrans, :), 'MarkerFaceColor', colorBehv(iTrans, :), 'LineStyle', 'none'); hold on; 
    plot(barPos_i(2), transAvg(2), 'Marker', '.', 'MarkerSize', 10, 'Color', colorBehv(iTrans, :), 'MarkerFaceColor', colorBehv(iTrans, :), 'LineStyle', 'none'); hold on; 
    errorbar(barPos_i(1), transAvg(1), transSem(1), 'Color', colorBehv(iTrans, :), 'Marker', '.', 'MarkerFaceColor', colorBehv(iTrans, :), 'LineStyle', '-', 'LineWidth', 2, 'MarkerSize', 15); 
    errorbar(barPos_i(2), transAvg(2), transSem(2), 'Color', colorBehv(iTrans, :), 'Marker', '.', 'MarkerFaceColor', colorBehv(iTrans, :), 'LineStyle', '-', 'LineWidth', 2, 'MarkerSize', 15); 

    [h, p, ci, stats] = ttest(data_i(:, 1), data_i(:, 2), 'Tail', tail_id)
    disp(['==============', num2str(iTrans), '==============']);
    disp(['t=', num2str(stats.tstat, '%4.3f'), ', p=', num2str(p, '%4.3f')])
end
set(gca, 'FontSize', 20, 'FontWeight', 'Bold', 'LineWidth', 2);
set(gca, 'XTick', '', 'XTickLabel', '');
xlim([0.3, 1.9]);
ylim([0.5, 2]);
box off;






