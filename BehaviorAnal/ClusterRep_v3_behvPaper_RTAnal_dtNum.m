% ClusterRep_v3_RTAnal_dtNum
% write by rxj @ 11/03/2020
% =========================================================================
% ClusterRep_v2_RTAnal_dtNum
% write by rxj @ 09/26/2020
% Control analysis: test the influence of distrator numbers on RTs

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
                 'sb_21_m_22'}; 
    subjLab   = {'hsp1', 'hr2', 'pxy3', 'wyx4', 'lml5', 'lf6', 'md7', 'srz8', 'fsq9', 'caq10', ...
                 'xjm11', 'cjs12', 'xxx13', 'wsq14', 'zzm15', 'lsy16', 'man17', 'zxy18', 'cyh19', 'szn20', ...
                 'sb21'};  
elseif isequal(ExpWord, 'cueRand')
    subj_list = {'zyh_1_f_21', 'why_2_m_18', 'cr_3_f_19', 'zyx_4_f_25', 'wym_5_m_23', 'wd_6_f_18', 'lyh_7_f_21', 'zr_8_f_24', 'zyh_9_m_19', 'zzy_10_f_19', ...
                 'smq_11_f_25', 'sz_12_f_20', 'lzy_13_m_22', 'yxy_14_f_18', 'zxl_15_f_18', 'wqh_16_m_23', 'zxj_17_m_21'}; % , 'xr_18_f_23'
    subjLab   = {'zyh1', 'why2', 'cr3', 'zyx4', 'wym5', 'wd6', 'lyh7', 'zr8', 'lyh9', 'zzy10', ...
                 'smq11', 'sz12', 'lzy13', 'yxy14', 'zxl15', 'wqh16', 'zxj17'}; % , 'xr18'
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
folder     = pwd;
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
RTtransTp_dtNum = zeros(subLen, 2, length(dtrNums));
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
    %%% within vs. between cluster transition RTs across different distractor numbers
    dtNum = arrayfun(@(x) length(find(dt_nodes(x, :)~=0)), (1 : nTrials)');
    dtNum = dtNum - 1;
    for iDt = 1 : length(dtrNums)
        RTtransTp_dtNum(SubIdx, 1, iDt) = mean(respRT_Col(find(dtNum == dtrNums(iDt) & transStyle == 1 & errorIdx_Col ~= 1 & fast_Col ~= 1)));
        RTtransTp_dtNum(SubIdx, 2, iDt) = mean(respRT_Col(find(dtNum == dtrNums(iDt) & transStyle ~= 1 & errorIdx_Col ~= 1 & fast_Col ~= 1)));
    end
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

%% (1) random + hamiltonian: within- and between-cluster RTs change with blocks
[RTtransTp_dtNumAvg, RTtransTp_dtNumSem] = Mean_and_Se(RTtransTp_dtNum, 1); 
RTtransTp_dtNumAvg = squeeze(RTtransTp_dtNumAvg);
RTtransTp_dtNumSem = squeeze(RTtransTp_dtNumSem);
figure('Position', [100 100 400 300]), clf;
errorbar(1 : 1 : length(dtrNums), RTtransTp_dtNumAvg(1, :), RTtransTp_dtNumSem(1, :), 'Color', colorSet(1, :), 'LineStyle', '-', 'LineWidth', 3); hold on;
errorbar(1 : 1 : length(dtrNums), RTtransTp_dtNumAvg(2, :), RTtransTp_dtNumSem(2, :), 'Color', colorSet(2, :), 'LineStyle', '-', 'LineWidth', 3); hold on;
set(gca, 'FontSize', 20, 'FontWeight', 'Bold', 'LineWidth', 2);
set(gca, 'XTick', 1 : 1 : length(dtrNums), 'XTickLabel', '');
xlim([0.5, length(dtrNums)+0.5]);
box off;
xlabel('Distractor Numbers', 'FontSize', 20);
ylabel('RTs (s)', 'FontSize', 20);

%%% repeated measurement ANOVA for gamma 
%  2 (cycle) by 2 (numerosity) repeated-measures ANOVAs
ID     = reshape(repmat((1 : 1 : subLen), length(dtrNums) * 2, 1), subLen * length(dtrNums) * 2, 1);
dtCol = repmat((1 : 1 : length(dtrNums))', subLen, 2); %% 7 blocks and 2 transitions (withih vs. between)
dtCol = reshape(dtCol', size(dtCol, 1)*size(dtCol, 2), 1);
transCol = repmat([1; 2], subLen * length(dtrNums), 1);
model1 = [1, 0, 0;... % main effect: subject
          0, 1, 0;... % main effect: distractor numbers
          0, 0, 1; ...% main effect: transitions
          0, 1, 1;...
          1, 1, 0; ...
          1, 0, 1];  
model = [1, 0, 0;... % main effect: subject
         0, 1, 0;... % main effect: distractor numbers
         0, 0, 1; ...% main effect: transitions
         0, 1, 1]; 
RTtest = reshape(RTtransTp_dtNum, [subLen, 2 * length(dtrNums)]);
RTtest = reshape(RTtest', [subLen * 2 * length(dtrNums), 1]);
[p_gam, tbl_gam, stats_gam] = anovan(RTtest, {ID, dtCol, transCol}, 'model', model, 'random', [1], 'sstype', 3, 'varnames', {'Subjects', 'distNo', 'TransType'}); % 'model' as 'interaction' is completely the same as model1



