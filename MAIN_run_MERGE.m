%% Load in Data
clear;
addpath(genpath('.'));

file = 'H:\Data\Other\MotionCorrected\Set_8_Shen_Shadlen\DNMF_Out_X6_from_full_v5.mat';

load(file);
options.SKEW_THR = 1.5;
options.SM_AMOUNT = 3;
options.DETREND_FRAMES = 900;
Cs2 = prep_FR(Cs, options.SM_AMOUNT, options.DETREND_FRAMES);
    
TEMPORAL_CORR_THR = options.temporalCorrThr;            % This is the big option to change to adjust how much merging happens
spatialOverlapThr = 0;
max_thr = options.maxthr;
min_sz = min(options.sizeRange);

skew = skewness(Cs2,[],2);
valid = skew>options.SKEW_THR;
Cs = Cs2(valid,:);
cROIs = cROIs(:,valid);
patchID = patchID(valid);

%% Merge similar ROIs
[Cs2, cROIs2] = patchMerge(Cs, cROIs, patchID, TEMPORAL_CORR_THR, spatialOverlapThr, options.SM_AMOUNT, options.DETREND_FRAMES, max_thr, min_sz, dimensions);


%% Visualize results
clf;
subplot(1,3,1);
imagescc(max(cROIs2,[],2));
axis square;
XL = xlim;
YL = ylim;

subplot(1,3,2);
hold on;
sz = sum(cROIs2>0);
[~,order] = sort(sz,'descend');
nROIs = size(cROIs2,2);
colors = distinguishable_colors(10);

for i=1:nROIs
    this = reshape(cROIs2(:,order(i)),dimensions);
%     maxx(i) = max(cROIs3(:,order2(i)))*max(Cs3(i,:));
    [bb,aa] = find(this>0);
    c = boundary(bb,aa,1);
    p = patch(aa(c),bb(c),colors(mod(i-1,10)+1,:));
    p.FaceAlpha = 0.75;
end
set(gca,'ydir','reverse');
xlim(XL);
ylim(YL);
axis square;

subplot(1,3,3);
imagesc(zscore(Cs2,[],2));
caxis([0 5]);
axis square;

%% Save to file
cROIs0 = cROIs;
cROIs = cROIs2;
Cs0 = Cs;
Cs = Cs2;
folder = fileparts(file);
save(fullfile(folder,'DNMF_Merged.mat'), 'cROIs0', 'Cs0', 'cROIs', 'Cs', 'options', 'dimensions', '-v7.3');
