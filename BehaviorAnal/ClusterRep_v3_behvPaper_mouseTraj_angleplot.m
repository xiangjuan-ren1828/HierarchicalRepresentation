% ClusterRep_v3_mouseTraj_angleplot
% write by rxj @ 07/08/2021
% angle change within a trial and across trials
%
% The angle data was calculated in jupyter lab and saved as concentricCir_pd_subj.csv

%%
clear
clc

%%
addpath('tight_subplot/');
addpath(genpath('HierarchicalCluster/'));
addpath('CircStat2012a/');

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
end
subjNum   = 1 : length(subj_list);
subj_listBv = subj_list;
subLen      = length(subj_list);

%% Experiment parameters
folder = pwd;
if isequal(ExpWord, 'cue')
    bhvDataDir = [folder, '/HierarchicalCluster/FormalExp-LynnNetwork-Results/'];
elseif isequal(ExpWord, 'explicit')
    bhvDataDir = [folder, '/HierarchicalCluster/FormalExp-ExplicitLearning-LynnNetwork-Results/'];
elseif isequal(ExpWord, 'cueRand')
    bhvDataDir = [folder, '/HierarchicalCluster/FormalExp-LynnNetwork-ImplicitRandom-Results/'];
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

% monitor parameters
scrWidthCM  = 47.9;
scrHeightCM = 26.9;
scrWidthPX  = 2048;
scrHeightPX = 1152;
winRect     = [0, 0, scrWidthPX, scrHeightPX];
visualAngle    = 9;  % 9 degree (diameter)
visualAngleCir = (visualAngle ./ 180) * pi;
visualDist     = 86; % 86 cm
dispDist       = 2 * visualDist .* tan(visualAngleCir / 2);
RChunk         = dispDist / 2;
imgWid         = 1/2 * RChunk;

cm2pxV      = scrHeightPX/scrHeightCM;
cm2pxH      = scrWidthPX/scrWidthCM;
[centerX, centerY] = RectCenter(winRect);
RChunkPX    = RChunk * cm2pxV;

imgSide = 1/2 * imgWid * cm2pxV; % or 100 pixel
imgSize = imgSide * 2;
angCir  = 0 : pi/50 : 2 * pi;
xCir    = RChunk * cm2pxH * cos(angCir) + centerX;
yCir    = RChunk * cm2pxV * sin(angCir) + centerY;    

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

%% read csv and plot 
concentricCir_pd_subj = readtable([bhvDataDir, 'concentricCir_pd_subj.csv']); %% add your own path

%%
subNo_col     = concentricCir_pd_subj.subNo;
trialCnt_col  = concentricCir_pd_subj.trialCnt;
blockNo_col   = concentricCir_pd_subj.blockNo;
trialType_col = concentricCir_pd_subj.trialType;
sampleCnt_col = concentricCir_pd_subj.sampleCnt;
radiusCon_col = concentricCir_pd_subj.radiusCon;
crossX_col    = concentricCir_pd_subj.crossX;
crossY_col    = concentricCir_pd_subj.crossY;
tgtAng_col    = concentricCir_pd_subj.tgtAng;
dtrAng1_col   = concentricCir_pd_subj.dtrAng1;
dtrAng2_col   = concentricCir_pd_subj.dtrAng2;
dtrAng3_col   = concentricCir_pd_subj.dtrAng3;
tgtXnorm_col  = concentricCir_pd_subj.tgtXnorm;
tgtYnorm_col  = concentricCir_pd_subj.tgtYnorm;
dtrXnorm1_col = concentricCir_pd_subj.dtrXnorm1;
dtrYnorm1_col = concentricCir_pd_subj.dtrYnorm1;
dtrXnorm2_col = concentricCir_pd_subj.dtrXnorm2;
dtrYnorm2_col = concentricCir_pd_subj.dtrYnorm2;
dtrXnorm3_col = concentricCir_pd_subj.dtrXnorm3;
dtrYnorm3_col = concentricCir_pd_subj.dtrYnorm3;
RT_col        = concentricCir_pd_subj.RT;
transType_col = concentricCir_pd_subj.transType;
dtrNoCnt_col  = concentricCir_pd_subj.dtrNoCnt; %% 1, 2, 3
trialMark2_col= concentricCir_pd_subj.trialMark2;
%%
circle_list = 0.02 : 0.02 : 0.98;
angleMatrix_subj = zeros(nTrials, length(circle_list), subLen);
itemPosX_subj    = zeros(nTrials, 4, subLen);
itemPosY_subj    = zeros(nTrials, 4, subLen);
crossXY_subj     = zeros(nTrials, length(circle_list), 2, subLen);
RT_subj    = zeros(nTrials, subLen);
trans_subj = zeros(nTrials, subLen);
dtrNo_subj = zeros(nTrials, subLen);
mark_subj  = zeros(nTrials, subLen);
for SubIdx = 1 : subLen
    subNo_i    = find(subNo_col == (SubIdx - 1));
    %%% angle
    tgtAng_i   = tgtAng_col(subNo_i);
    tgtAng_mat = reshape(tgtAng_i, length(circle_list), nTrials);
    angleMatrix_subj(:, :, SubIdx) = tgtAng_mat'; %%% trial * radius
    %%% pixel positions of target and distractors
    tgtXnorm_i  = tgtXnorm_col(subNo_i);
    tgtYnorm_i  = tgtYnorm_col(subNo_i);
    dtrXnorm1_i = dtrXnorm1_col(subNo_i);
    dtrYnorm1_i = dtrYnorm1_col(subNo_i);
    dtrXnorm2_i = dtrXnorm2_col(subNo_i);
    dtrYnorm2_i = dtrYnorm2_col(subNo_i);
    dtrXnorm3_i = dtrXnorm3_col(subNo_i);
    dtrYnorm3_i = dtrYnorm3_col(subNo_i);
    tgtXnorm_t  = reshape(tgtXnorm_i, length(circle_list), nTrials);
    tgtYnorm_t  = reshape(tgtYnorm_i, length(circle_list), nTrials);
    dtrXnorm1_t = reshape(dtrXnorm1_i, length(circle_list), nTrials);
    dtrYnorm1_t = reshape(dtrYnorm1_i, length(circle_list), nTrials);
    dtrXnorm2_t = reshape(dtrXnorm2_i, length(circle_list), nTrials);
    dtrYnorm2_t = reshape(dtrYnorm2_i, length(circle_list), nTrials);
    dtrXnorm3_t = reshape(dtrXnorm3_i, length(circle_list), nTrials);
    dtrYnorm3_t = reshape(dtrYnorm3_i, length(circle_list), nTrials);
    itemPosX_subj(:, :, SubIdx) = [tgtXnorm_t(1, :); dtrXnorm1_t(1, :); dtrXnorm2_t(1, :); dtrXnorm3_t(1, :)]';
    itemPosY_subj(:, :, SubIdx) = [tgtYnorm_t(1, :); dtrYnorm1_t(1, :); dtrYnorm2_t(1, :); dtrYnorm3_t(1, :)]';
    %%% interpolation X and Y of the mouse trajectory
    crossX_i = crossX_col(subNo_i);
    crossY_i = crossY_col(subNo_i);
    crossX_t = reshape(crossX_i, length(circle_list), nTrials);
    crossY_t = reshape(crossY_i, length(circle_list), nTrials);
    crossXY_subj(:, :, 1, SubIdx) = crossX_t';
    crossXY_subj(:, :, 2, SubIdx) = crossY_t';
    %%% RT
    RT_i     = RT_col(subNo_i);
    RT_trial = reshape(RT_i, length(circle_list), nTrials);
    RT_subj(:, SubIdx) = (RT_trial(1, :))';
    %%% transition
    transType_i = transType_col(subNo_i);
    transType_t = reshape(transType_i, length(circle_list), nTrials);
    trans_subj(:, SubIdx) = (transType_t(1, :))';
    %%% distractor number
    dtrNoCnt_i  = dtrNoCnt_col(subNo_i);
    dtrNoCnt_t  = reshape(dtrNoCnt_i, length(circle_list), nTrials);
    dtrNo_subj(:, SubIdx) = (dtrNoCnt_t(1, :))';
    %%% trial mark
    trialMark2_i = trialMark2_col(subNo_i);
    trialMark2_t = reshape(trialMark2_i, length(circle_list), nTrials);
    mark_subj(:, SubIdx) = (trialMark2_t(1, :))';
end

%% plot all trials for each condition
for SubIdx = 1% : subLen
    subTit     = subjLab{SubIdx};
    angle_iSub = angleMatrix_subj(:, :, SubIdx); %% trial * radius
    RT_iSub    = RT_subj(:, SubIdx);
    trans_iSub = trans_subj(:, SubIdx);
    mark_iSub  = mark_subj(:, SubIdx);
    
    %%% 4 heatmap: (1) correct before cue; (2) correct after cue; (3) incorrect before cue; (4) incorrect after cue
    figure('Position', [100 100 1000 220]), clf;
    ha = tight_subplot(1, 4, [.05 .05], [.1 .02], [.05 .05]);
    ha_i = 1;
    for i = 1 : 4
        mark_i  = find(mark_iSub == (i - 1));
        angle_i = angle_iSub(mark_i, :);
        if ~isempty(mark_i)
            iT_del = [];
            for iT = 1 : size(angle_i, 1)
                angle_i_iT = angle_i(iT, :);
                if any(isnan(angle_i_iT))
                    iT_del = [iT_del; iT];
                end
            end
            mark_i(iT_del)     = [];
            angle_i(iT_del, :) = [];
        end  
        trans_i = trans_iSub(mark_i, :);
        trans_w = find(trans_i == 0 | trans_i == 1 | trans_i == 2); %%% within-cluster transition
        trans_b = find(trans_i == 3); %%% between-cluster transition
        
        axes(ha(ha_i));
        yyaxis left
        imagesc(circle_list, [1, length(mark_i)], angle_i); hold on;
        %[colorCohr]=cbrewer('div', 'RdBu', 100, 'cubic');
        [colorCohr]=cbrewer('div', 'PuOr', 100, 'cubic'); %'PRGn'
        colormap(colorCohr(end : -1 : 1, :)); % colormap(fliplr(colorCohr));
        colorbar;
        %caxis([-10 25]);
        axis xy;
        set(gca, 'LineWidth', 2);
        set(gca, 'FontSize', 15, 'FontWeight', 'bold', 'FontName', 'Arial');
        box off;
        
        yyaxis right
        ylim([-4, 4]);
        if ~isempty(trans_w)
            plot(circle_list, circ_mean(angle_i(trans_w, :), [], 1), 'Color', 'r', 'LineStyle', '-', 'LineWidth', 3); hold on;
        end
        if ~isempty(trans_b)
            plot(circle_list, circ_mean(angle_i(trans_b, :), [], 1), 'Color', 'b', 'LineStyle', '-', 'LineWidth', 3); hold on;
        end
        %title(subTit, 'FontSize', 14);
        ha_i = ha_i + 1;
    end
end

%% color settings
colorSet = [249, 183, 176; ...
            248, 218, 172; ...
            84, 185, 211; ...
            184, 204, 225; ...
            138, 170, 51; ...
            210, 234, 200; ...
            198, 127, 192; ...
            219, 204, 226] ./ [255, 255, 255];
        
%% plot the angle mean and variance for individual subject and for each distractor number, merging all 1500 trials
stat_subj = zeros(4, length(circle_list), 3, 2, subLen); 
for SubIdx = 1 : subLen
    subTit     = subjLab{SubIdx};
    angle_iSub = angleMatrix_subj(:, :, SubIdx);
    mark_iSub  = mark_subj(:, SubIdx);
    dtrNo_iSub = dtrNo_subj(:, SubIdx);
    %%% data re-organization
    stat_iSub = nan(4, length(circle_list), 3, 2); %% 2 means 'mean' and 'variance'
    for j = 1 : 3      %% 3 distractor number
        for i = 1 : 4  %% 4 trial types: (1) correct before cue; (2) correct after cue; (3) incorrect before cue; (4) incorrect after cue
            mark_i  = find((mark_iSub == (i - 1)) & (dtrNo_iSub == j));
            angle_i = angle_iSub(mark_i, :);
            %%% check whether there is NAN within a trial
            if ~isempty(mark_i)
                iT_del = [];
                for iT = 1 : size(angle_i, 1)
                    angle_i_iT = angle_i(iT, :);
                    if any(isnan(angle_i_iT))
                        iT_del = [iT_del; iT];
                    end
                end
                mark_i(iT_del)     = [];
                angle_i(iT_del, :) = [];
            end
            if ~isempty(mark_i)
                stat_iSub(i, :, j, 1) = circ_mean(angle_i, [], 1);
                stat_iSub(i, :, j, 2) = circ_var(angle_i, [], [], 1);
            end
            
        end
    end
    stat_subj(:, :, :, :, SubIdx) = stat_iSub;
    
    figure('Position', [100 100 800 450]), clf;
    ha = tight_subplot(2, 3, [.05 .05], [.1 .02], [.05 .05]);
    ha_i = 1;
    for m = 1 : 2          %% mean (1st row) and variance (2nd row)
        for j = 1 : 3      %% 3 distractor number
            axes(ha(ha_i));
            for i = 1 : 4  %% 4 trial types: (1) correct before cue; (2) correct after cue; (3) incorrect before cue; (4) incorrect after cue
                circ_stat = stat_iSub(i, :, j, m);
                if ~all(isnan(circ_stat))
                    plot(circle_list, circ_stat, 'Color', colorSet(i, :), 'LineStyle', '-', 'LineWidth', 3); hold on;
                end
            end
            if m == 1
                ylim([-4, 4]); %% mean
            elseif m == 2
                ylim([0, 1]);  %% variance
            end
            set(gca, 'LineWidth', 2);
            set(gca, 'FontSize', 15, 'FontWeight', 'bold', 'FontName', 'Arial');
            box off;
            ha_i = ha_i + 1;
        end
    end
end
%% plot the mouse trajectory for individual subject; corresponding to the above figure
% rotate the target to (1, 0);s
%itemPosX_subj    = zeros(nTrials, 4, subLen);
%itemPosY_subj    = zeros(nTrials, 4, subLen);
data_id = input('which data to plot, 0-raw, 1-normalized: ');
SubIdx = 1;
subjBv  = subj_listBv{SubIdx};
subTit  = subjLab{SubIdx};
subjDir = [bhvDataDir, subjBv, '/'];
mouseTraj_all = [];
for iBlock = 1 : nBlock
    load([subjDir, subjBv, 'mouseTraj_trials_in_Blc', num2str(iBlock), '.mat'], 'mouseTraj_trials');
    mouseTraj_all = [mouseTraj_all; mouseTraj_trials];
end
load([subjDir, subjBv, 'clusterResult_Blocks.mat'], 'clusterResult');
itemXraw = clusterResult(:, imgPosX); % imgPosX  = 17 : 20;
itemYraw = clusterResult(:, imgPosY); % imgPosY  = 21 : 24;

zeroPos = [1, 0]; %% the x- and y-pos of 0 angle image
j = 1; %% only one distractor
angle_iSub = angleMatrix_subj(:, :, SubIdx);
mark_iSub  = mark_subj(:, SubIdx);
dtrNo_iSub = dtrNo_subj(:, SubIdx);
itemX_iSub = itemPosX_subj(:, :, SubIdx);
itemY_iSub = itemPosY_subj(:, :, SubIdx);
crossXY_iSub = crossXY_subj(:, :, :, SubIdx);
for i = 1 : 4  %% 4 trial types: (1) correct before cue; (2) correct after cue; (3) incorrect before cue; (4) incorrect after cue
    mark_i     = find((mark_iSub == (i - 1)) & (dtrNo_iSub == j));
    itemX_i    = itemX_iSub(mark_i, :);
    itemY_i    = itemY_iSub(mark_i, :);
    itemXraw_i = itemXraw(mark_i, :);
    itemYraw_i = itemYraw(mark_i, :);
    %%% mouse trajectory
    mouseTraj_i = mouseTraj_all(mark_i, :);
    crossXY_i   = crossXY_iSub(mark_i, :, :); %% length * radius * 2 (XY)
    
    if  i == 1 || i == 2 || i == 4
        figure('Position', [100 100 1500 700]), clf;
        ha = tight_subplot(3, 6, [.05 .05], [.1 .02], [.05 .05]);
    elseif i == 3
        figure('Position', [100 100 1500 800]), clf;
        ha = tight_subplot(4, 6, [.05 .05], [.1 .02], [.05 .05]);
    end
    if i == 1 || i == 2
        plot_len = 18;
    elseif i == 3 || i == 4
        plot_len = length(mark_i);
    end
    ha_i = 1;
    for iPlot = 1 : plot_len
        axes(ha(ha_i));
        if data_id == 0
            colorTmp = cbrewer('div', 'RdBu', length(mouseTraj_i{iPlot}(:, 1)), 'cubic'); % parula(length(mouseTraj_i{iPlot}(:, 1)));
            colorTmp = colorTmp(end : -1 : 1, :);
            plot(xCir, yCir, 'Marker', 'none', 'Color', [0.4, 0.4, 0.4], 'LineStyle', ':', 'LineWidth', 1); hold on;
            for iDot = 1 : length(mouseTraj_i{iPlot}(:, 1))
                plot(mouseTraj_i{iPlot}(iDot, 1), mouseTraj_i{iPlot}(iDot, 2), 'Color', colorTmp(iDot, :), 'Marker', '.', 'MarkerSize', 10); hold on;
            end
            plot(itemXraw_i(iPlot, 1), itemYraw_i(iPlot, 1), 'Marker', '+', 'MarkerSize', 15, 'Color', 'r', 'LineWidth', 2); hold on;
            plot(itemXraw_i(iPlot, 2), itemYraw_i(iPlot, 2), 'Marker', '*', 'MarkerSize', 15, 'Color', 'b', 'LineWidth', 2); hold on;
            
        elseif data_id == 1
            %%% -----------------------------------------------------------
            %%% ratote target to (1,0); all of the other distractors and
            %%% mouse landing positions were rotated according to the angle
            %%% between target and (1, 0)
            %%% -----------------------------------------------------------
            %%% calculate the angle between target and zero angle
            Col1 = ([itemX_i(iPlot, 1), itemY_i(iPlot, 1)])' - [0; 0];
            Col2 = zeroPos' - [0; 0];
            cosTheta = dot(Col1, Col2) / (norm(Col1) * norm(Col2)); %% for Col1 and Col2, each column denotes a unique point
            % Theta = acos(cosTheta) * 180 / pi; % unit in degree
            if itemY_i(iPlot, 1) <= 0
                thetaTgt = acos(cosTheta); %% if the target's Y-axis is smaller than the zero pionts, then the rotation should be counterclockwise, that is the angle should be positive
            else
                thetaTgt = -acos(cosTheta); %% otherwise, the rotation should be clockwise
            end
            %%% all of the data points multiplied by the rotation
            %%% matrix
            rotMat = [cos(thetaTgt), -sin(thetaTgt); ...
                      sin(thetaTgt), cos(thetaTgt)];
            posXY = squeeze(crossXY_i(iPlot, :, :));   
            posXY_d = posXY' - [0; 0];
            rotXY   = rotMat * posXY_d + [0; 0];
            tgtPos_rot = rotMat * ([itemX_i(iPlot, 1), itemY_i(iPlot, 1)])';
            dtrPos_rot = rotMat * ([itemX_i(iPlot, 2), itemY_i(iPlot, 2)])';

            %%% plotting
            colorTmp = cbrewer('div', 'RdBu', length(circle_list), 'cubic'); % parula(length(mouseTraj_i{iPlot}(:, 1)));
            colorTmp = colorTmp(end : -1 : 1, :);
            plot(cos(angCir), sin(angCir), 'Marker', 'none', 'Color', [0.4, 0.4, 0.4], 'LineStyle', ':', 'LineWidth', 1); hold on;
            for iDot = 1 : length(circle_list)
                %plot(squeeze(crossXY_i(iPlot, iDot, 1)), squeeze(crossXY_i(iPlot, iDot, 2)), 'Color', colorTmp(iDot, :), 'Marker', '.', 'MarkerSize', 10); hold on;
                plot(rotXY(1, iDot), rotXY(2, iDot), 'Color', colorTmp(iDot, :), 'Marker', '.', 'MarkerSize', 10); hold on;
            end
            %plot(itemX_i(iPlot, 1), itemY_i(iPlot, 1), 'Marker', '+', 'MarkerSize', 15, 'Color', 'r', 'LineWidth', 2); hold on;
            %plot(itemX_i(iPlot, 2), itemY_i(iPlot, 2), 'Marker', '*', 'MarkerSize', 15, 'Color', 'b', 'LineWidth', 2); hold on;
            plot(tgtPos_rot(1), tgtPos_rot(2), 'Marker', '+', 'MarkerSize', 15, 'Color', 'r', 'LineWidth', 2); hold on;
            plot(dtrPos_rot(1), dtrPos_rot(2), 'Marker', '*', 'MarkerSize', 15, 'Color', 'b', 'LineWidth', 2); hold on;
        end
        set(gca, 'LineWidth', 2);
        set(gca, 'FontSize', 15, 'FontWeight', 'bold', 'FontName', 'Arial');
        axis square
        if data_id == 0
            xlim([600, 1400]);
            ylim([200, 1000]);
        elseif data_id == 1
            xlim([-1.4, 1.4]);
            ylim([-1.4, 1.4]);
        end
        box off;
        ha_i = ha_i + 1;
        
    end
end
%%
figure('Position', [50, 50, 190, 400]), clf;
%colormap(jet);
[colorCohr]=cbrewer('div', 'RdBu', 100, 'cubic');
colormap(colorCohr(end : -1 : 1, :)); % colormap(fliplr(colorCohr)); 
ax = axes;
% the following properties of colorbar was described with details at
% 'ColorBar Properties' (https://www.mathworks.com/help/matlab/ref/matlab.graphics.illustration.colorbar-properties.html)
c = colorbar(ax, 'Location', 'east',  'Ticks', '', 'TickLabels', '', 'TickLength', 0.02, 'FontSize', 28, 'FontWeight', 'Bold', 'Position', [0.46, 0.1, 0.1, 0.8]); %{'low', 'medium', 'high'}'FontAngle', 'italic',
%c.Label.String = 'Decoding Performance';
c.LineWidth = 1.5;
c.Label.FontSize = 28;
c.Label.FontWeight = 'Bold';
%c.Label.FontAngle  = 'italic';
ax.Visible = 'off';
%cfg.zlim   = [-0.02, 0.04];

%% angular distribution across trials at each radius, separately for 4 kinds of trials (1. correct before cue; 2. correct after cue; 3. incorrect before cue; 4. incorrect after cue)
% added by rxj @ 07/11/2021
for SubIdx = 1% : subLen
    subTit     = subjLab{SubIdx};
    angle_iSub = angleMatrix_subj(:, :, SubIdx);
    mark_iSub  = mark_subj(:, SubIdx);
    dtrNo_iSub = dtrNo_subj(:, SubIdx);
    %%% data re-organization
    stat_iSub = nan(4, length(circle_list), 3, 2); %% 2 means 'mean' and 'variance'
    for j = 3% : 3      %% 3 distractor number
        for i = 1 : 4  %% 4 trial types: (1) correct before cue; (2) correct after cue; (3) incorrect before cue; (4) incorrect after cue
            mark_i  = find((mark_iSub == (i - 1)) & (dtrNo_iSub == j));
            angle_i = angle_iSub(mark_i, :);
            %%% check whether there is NAN within a trial
            if ~isempty(mark_i)
                iT_del = [];
                for iT = 1 : size(angle_i, 1)
                    angle_i_iT = angle_i(iT, :);
                    if any(isnan(angle_i_iT))
                        iT_del = [iT_del; iT];
                    end
                end
                mark_i(iT_del)     = [];
                angle_i(iT_del, :) = [];
            end
            
            figure('Position', [100 100 1200 400]), clf;
            ha = tight_subplot(5, 10, [.05 .01], [.02 .05], [.05 .05]);
            r = 1.2;
            ha_i = 1;
            for iCir = 1 : length(circle_list)
                axes(ha(ha_i));
                x = 0;
                y = 0;
                viscircles([x, y], r, 'Color', colorSet(i, :), 'LineStyle', '-'); hold on;
                daspect([1,1,1]);
                
                PhaseDiffTmp = angle_i(:, iCir); % one column
                LineXsub     = x + r * cos(PhaseDiffTmp);
                LineYsub     = y + r * sin(PhaseDiffTmp);
                
                plot([x - r, x + r], [y, y], 'Color', 'k', 'LineStyle', '-', 'LineWidth', 0.5); hold on;
                plot([x, x], [y - r, y + r], 'Color', 'k', 'LineStyle', '-', 'LineWidth', 0.5); hold on;
                % single subject
                for iT = 1 : length(PhaseDiffTmp)
                    plot([x, LineXsub(iT, :)], [y, LineYsub(iT, :)], 'Color', [0.5, 0.5, 0.5], 'LineStyle', '-', 'LineWidth', 1); hold on;
                end
                xlim([-1.5, 1.5]);
                ylim([-1.5, 1.5]);
                
                meanTmp      = circ_mean(PhaseDiffTmp);
                LineX        = x + r * cos(meanTmp); % x and y of phase angle
                LineY        = y + r * sin(meanTmp);
                pCir         = circ_rtest(PhaseDiffTmp);
                % subject average
                if pCir < 0.05
                    plot([x, LineX], [y, LineY], 'Color', 'r', 'LineStyle', '-', 'LineWidth', 3); hold on;
                else
                    plot([x, LineX], [y, LineY], 'Color', 'k', 'LineStyle', '-', 'LineWidth', 3); hold on;
                end
                box off;
                title(['r=', num2str(circle_list(iCir), '%1.2f')])
                ha_i = ha_i + 1;
            end
            %set(gca, 'FontSize', 16);
            %set(gca, 'LineWidth', 2);
        end
    end
end


%% subjects average
stat_circmean_subjAvg = zeros(4, length(circle_list), 3);
for i = 1 : 4
    for iC = 1 : length(circle_list)
        for j = 1 : 3
            circ_mean_i = squeeze(stat_subj(i, iC, j, 1, :)); %stat_subj = zeros(4, length(circle_list), 3, 2, subLen);
            circ_mean_i(isnan(circ_mean_i)) = [];
            stat_circmean_subjAvg(i, iC, j) = circ_mean(circ_mean_i);
        end
    end
end
%stat_circmean_subjAvg = circ_mean(squeeze(stat_subj(:, :, :, 1, :)), [], 4);
stat_circvar_subjAvg  = nanmean(squeeze(stat_subj(:, :, :, 2, :)), 4);
stat_avg = zeros(4, length(circle_list), 3, 2); 
stat_avg(:, :, :, 1) = stat_circmean_subjAvg;
stat_avg(:, :, :, 2) = stat_circvar_subjAvg;

figure('Position', [100 100 800 450]), clf;
ha = tight_subplot(2, 3, [.05 .05], [.1 .02], [.05 .05]);
ha_i = 1;
for m = 1 : 2          %% mean (1st row) and variance (2nd row)
    for j = 1 : 3      %% 3 distractor number
        axes(ha(ha_i));
        for i = 1 : 4  %% 4 trial types: (1) correct before cue; (2) correct after cue; (3) incorrect before cue; (4) incorrect after cue
            circ_stat = stat_avg(i, :, j, m);
            plot(circle_list, circ_stat, 'Color', colorSet(i, :), 'LineStyle', '-', 'LineWidth', 3); hold on;
        end
%         if m == 1
%             ylim([-4, 4]);
%         elseif m == 2
%             ylim([0, 1]);
%         end
        set(gca, 'LineWidth', 2);
        set(gca, 'FontSize', 15, 'FontWeight', 'bold', 'FontName', 'Arial');
        box off;
        ha_i = ha_i + 1;
    end
end

%% angular distribution across subjects at each radius
% added by rxj @ 07/11/2021
for j = 3% : 3      %% 3 distractor number
    for i = 1 : 4  %% 4 trial types: (1) correct before cue; (2) correct after cue; (3) incorrect before cue; (4) incorrect after cue
        figure('Position', [100 100 1200 400]), clf;
        ha = tight_subplot(5, 10, [.05 .01], [.02 .05], [.05 .05]);
        r = 1.2;
        ha_i = 1;
        for iCir = 1 : length(circle_list)
            axes(ha(ha_i));
            x = 0;
            y = 0;
            viscircles([x, y], r, 'Color', colorSet(i, :), 'LineStyle', '-'); hold on;
            daspect([1,1,1]);
            
            PhaseDiffTmp = squeeze(stat_subj(i, iCir, j, 1, :)); % one column
            LineXsub     = x + r * cos(PhaseDiffTmp);
            LineYsub     = y + r * sin(PhaseDiffTmp);
            
            plot([x - r, x + r], [y, y], 'Color', 'k', 'LineStyle', '-', 'LineWidth', 0.5); hold on;
            plot([x, x], [y - r, y + r], 'Color', 'k', 'LineStyle', '-', 'LineWidth', 0.5); hold on;
            % single subject
            for iT = 1 : length(PhaseDiffTmp)
                plot([x, LineXsub(iT, :)], [y, LineYsub(iT, :)], 'Color', [0.5, 0.5, 0.5], 'LineStyle', '-', 'LineWidth', 1); hold on;
            end
            xlim([-1.5, 1.5]);
            ylim([-1.5, 1.5]);
            
            PhaseDiffTmp(isnan(PhaseDiffTmp)) = [];
            meanTmp      = circ_mean(PhaseDiffTmp);
            LineX        = x + r * cos(meanTmp); % x and y of phase angle
            LineY        = y + r * sin(meanTmp);
            pCir         = circ_rtest(PhaseDiffTmp);
            % subject average
            if pCir < 0.05
                plot([x, LineX], [y, LineY], 'Color', 'r', 'LineStyle', '-', 'LineWidth', 3); hold on;
            else
                plot([x, LineX], [y, LineY], 'Color', 'k', 'LineStyle', '-', 'LineWidth', 3); hold on;
            end
            box off;
            title(['r=', num2str(circle_list(iCir), '%1.2f')])
            ha_i = ha_i + 1;
        end
        %set(gca, 'FontSize', 16);
        %set(gca, 'LineWidth', 2);
    end
end



%% plot the angle mean and variance for individual subject, seperately for learning and testing session
%%% 1st row: first 700 trials (learning session); 2nd row: 800 trials (testing)
for SubIdx = 1%1 : subLen
    subTit     = subjLab{SubIdx};

    %%% 4 heatmap: (1) correct before cue; (2) correct after cue; (3) incorrect before cue; (4) incorrect after cue
    figure('Position', [100 100 1000 450]), clf;
    ha = tight_subplot(2, 4, [.05 .05], [.1 .02], [.05 .05]);
    ha_i = 1;
    for j = 1 : 2
        if j == 1     %% learning
            angle_iSub = angleMatrix_subj(1 : rndTrial, :, SubIdx);
            RT_iSub    = RT_subj(1 : rndTrial, SubIdx);
            trans_iSub = trans_subj(1 : rndTrial, SubIdx);
            mark_iSub  = mark_subj(1 : rndTrial, SubIdx);
            
        elseif j == 2 %% testing
            angle_iSub = angleMatrix_subj((rndTrial + 1) : end, :, SubIdx);
            RT_iSub    = RT_subj((rndTrial + 1) : end, SubIdx);
            trans_iSub = trans_subj((rndTrial + 1) : end, SubIdx);
            mark_iSub  = mark_subj((rndTrial + 1) : end, SubIdx);
            
        end
        for i = 1 : 4
            mark_i  = find(mark_iSub == (i - 1));
            angle_i = angle_iSub(mark_i, :);
            if ~isempty(mark_i)
                iT_del = [];
                for iT = 1 : size(angle_i, 1)
                    angle_i_iT = angle_i(iT, :);
                    if any(isnan(angle_i_iT))
                        iT_del = [iT_del; iT];
                    end
                end
                mark_i(iT_del)     = [];
                angle_i(iT_del, :) = [];
            end
            trans_i = trans_iSub(mark_i, :);
            trans_w = find(trans_i == 0 | trans_i == 1 | trans_i == 2); %%% within-cluster transition
            trans_b = find(trans_i == 3); %%% between-cluster transition
            
            axes(ha(ha_i));
            yyaxis left
            ylim([-4, 4]);
            if ~isempty(trans_w)
                plot(circle_list, circ_mean(angle_i(trans_w, :), [], 1), 'Color', colorSet(1, :), 'LineStyle', '-', 'LineWidth', 3); hold on;
            end
            if ~isempty(trans_b)
                plot(circle_list, circ_mean(angle_i(trans_b, :), [], 1), 'Color', colorSet(1, :), 'LineStyle', ':', 'LineWidth', 3); hold on;
            end
            ylim([-4, 4]);
            ax = gca;
            ax.YColor = colorSet(1, :);
        
            yyaxis right
            ylim([0, 1]);
            if ~isempty(trans_w)
                plot(circle_list, circ_var(angle_i(trans_w, :), [], [], 1), 'Color', colorSet(3, :), 'LineStyle', '-', 'LineWidth', 3); hold on;
            end
            if ~isempty(trans_b)
                plot(circle_list, circ_var(angle_i(trans_b, :), [], [], 1), 'Color', colorSet(3, :), 'LineStyle', ':', 'LineWidth', 3); hold on;
            end
            ylim([0, 1]);
            ax = gca;
            ax.YColor = colorSet(3, :);
            
            set(gca, 'LineWidth', 2);
            set(gca, 'FontSize', 15, 'FontWeight', 'bold', 'FontName', 'Arial');
            box off;
            ha_i = ha_i + 1;
            
        end
    end
end































































