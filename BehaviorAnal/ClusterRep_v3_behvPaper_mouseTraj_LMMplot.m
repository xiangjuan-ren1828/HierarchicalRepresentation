% ClusterRep_v3_mouseTraj_LMMplot
% write by rxj @ 07/14/2021
% plot the linear mixed-effect model or GLM results

%%
clear
clc

%%
addpath('tight_subplot/');
addpath(genpath('HierarchicalCluster/'));

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
folder = '/Users/cdlab/Nextcloud/MATLAB/MyExperiment/HierarchicalCluster'; %pwd;
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
colorSets_d = colorSets;
colorSets_d(2, :) = [];
%circle_list = 0.02 : 0.02 : 0.98;
circle_list = 0 : 1/60 : 1.5;

%% define the dependent variables and data kind
regressorWord_id  = 2;
regressorWord_all = {'TrlCnt-TrlType', 'ClsId-TransType'}; % the first 2 parameters
regressorWord     = regressorWord_all{regressorWord_id};

dataFrom_id  = 3;
dataFrom_all = {'allTrials', 'learning', 'testing'}; % 'allTrials_learning', 'allTrials_testing',
dataFrom     = dataFrom_all{dataFrom_id};

depVar_id  = 1;
depVar_all = {'absAng', 'cumsumDist', 'insSpd'};
depVar     = depVar_all{depVar_id};

if isequal(dataFrom, 'allTrials')
    nParams = 10;
    colorUsed = colorSets;
else
    if isequal(regressorWord, 'TrlCnt-TrlType')
        nParams = 9;
        colorUsed = colorSets_d;
    elseif isequal(regressorWord, 'ClsId-TransType')
        nParams = 10;
        colorUsed = colorSets;
    end
end
if isequal(depVar, 'absAng')
    markPos = 0.49 : -0.01 : 0.40;
elseif isequal(depVar, 'cumsumDist')
    markPos = 0.6 : -0.01 : 0.51;
end

%% Analysis 
lme_coef_exp = zeros(nParams, length(circle_list), 3);
lme_pval_exp = zeros(nParams, length(circle_list), 3);
for i = 1 : 3
    if i == 1
        bhvDataDir = [folder, '/FormalExp-LynnNetwork-Results/angleCal/'];
    elseif i == 2
        bhvDataDir = [folder, '/FormalExp-ExplicitLearning-LynnNetwork-Results/angleCal/'];
    elseif i == 3
        bhvDataDir = [folder, '/FormalExp-LynnNetwork-ImplicitRandom-Results/angleCal/'];
    end
    %%% read data
    if isequal(regressorWord, 'TrlCnt-TrlType')
        lme_coef = readtable([bhvDataDir, 'lme_coef_', dataFrom, '_', depVar, '.csv']); %% add your own path
        lme_pval = readtable([bhvDataDir, 'lme_pval_', dataFrom, '_', depVar, '.csv']);
    elseif isequal(regressorWord, 'ClsId-TransType')
        lme_coef = readtable([bhvDataDir, 'lme_coef_clsIDtransType_', dataFrom, '_', depVar, '.csv']); %% add your own path
        lme_pval = readtable([bhvDataDir, 'lme_pval_clsIDtransType_', dataFrom, '_', depVar, '.csv']);
    end
    
    lme_coef(:, 1) = []; % paramsLen * 49 (concentric circle)
    lme_coef = table2array(lme_coef);
    lme_pval(:, 1) = []; 
    lme_pval = table2array(lme_pval);
    
    lme_coef_exp(:, :, i) = lme_coef;
    lme_pval_exp(:, :, i) = lme_pval;
end

%% plotting 
figure('Position', [100 100 1000 300]), clf;
ha = tight_subplot(1, 3, [.05 .05], [.1 .02], [.05 .05]);
ha_i = 1;
for i = 1 : 3
    coef_i = lme_coef_exp(:, :, i);
    pval_i = lme_pval_exp(:, :, i);

    axes(ha(ha_i));
    %%
    for j = 1 : size(coef_i, 1)
        plot(circle_list, coef_i(j, :), 'Color', colorUsed(j, :), 'LineStyle', '-', 'LineWidth', 2); hold on;
        [~, ~, ~, adj_p] = fdr_bh((pval_i(j, :))', 0.05, 'pdep'); % method: 'dep', 'pdep'
        
        for iCir = 1 : length(circle_list)
            pval_j = adj_p(iCir);
            if pval_j < 0.01
                plot(circle_list(iCir), markPos(j), 'Marker', '.', 'MarkerSize', 8, 'Color', colorUsed(j, :), 'MarkerFaceColor', colorUsed(j, :), 'LineStyle', 'none', 'LineWidth', 1.5); hold on;
            elseif pval_j < 0.05 && pval_j >= 0.01
                plot(circle_list(iCir), markPos(j), 'Marker', '>', 'MarkerSize', 4, 'Color', colorUsed(j, :), 'MarkerFaceColor', colorUsed(j, :), 'LineStyle', 'none', 'LineWidth', 1.5); hold on;
            end
            
        end
    end
    %xlim([0, 1]);
    plot(xlim, [0, 0], 'k--', 'LineWidth', 1); hold on;
    if isequal(regressorWord, 'TrlCnt-TrlType')
        if isequal(depVar, 'absAng')
            ylim([-0.3, 0.5]);
        elseif isequal(depVar, 'cumsumDist')
            ylim([-1, 0.6]);
        end
    end
    axis xy;
    set(gca, 'LineWidth', 2);
    set(gca, 'FontSize', 15, 'FontWeight', 'bold', 'FontName', 'Arial');
    box off;
    ha_i = ha_i + 1;
end

















