function [Cs, cROIs] = patchMerge(Cs0, cROIs0, patchID, temporalCorrThr, spatialOverlapThr, SM_AMOUNT, DETREND_FRAMES, max_pct, min_sz, fullSize)

    if(nargin<4)
        temporalCorrThr = 0.9;
    end
    if(nargin<5)
        spatialOverlapThr = 0;
    end
    if(nargin<6)
        SM_AMOUNT = 3;
    end
    if(nargin<7)
        DETREND_FRAMES = 900;
    end
    if(nargin<8)
        max_pct = 0.2;
    end
    if(nargin<9)
        min_sz = 20;
    end
    if(nargin<10)
        fullSize = [512 512];
    end
    
    
    CROIS = cell(1,max(patchID));
    CS = cell(max(patchID),1);
    % Loop over every patch and merge things within there. No normalization
    normalize_option = false;
    activity_option = false;
    
    
    fprintf('\nMerging within patches...\n');
    for i_patch = 1:max(patchID)
        these = (patchID == i_patch);
        if(sum(these)==0)
            continue;
        end
        theseROIs = cROIs0(:,these);
        theseCs = Cs0(these,:);
        
        [cROIs2, Cs2] = mergeROIs3(theseCs, theseROIs, temporalCorrThr, spatialOverlapThr, fullSize, normalize_option, activity_option, SM_AMOUNT, DETREND_FRAMES);
        [cROIs2b, Cs2b] = mergeROIs3(Cs2, cROIs2, temporalCorrThr, spatialOverlapThr, fullSize, normalize_option, activity_option, SM_AMOUNT, DETREND_FRAMES);

        CROIS{i_patch} = cROIs2b;
        CS{i_patch} = Cs2b;
    end

    Cs2 = cell2mat(CS);
    cROIs2 = cell2mat(CROIS);

    %% Now merge all ROIs with normalization
    normalize_option = true;
    fprintf('Global merging...\n');
    [cROIs3, Cs3] = mergeROIs3(Cs2, cROIs2, temporalCorrThr, spatialOverlapThr, fullSize, normalize_option, activity_option, SM_AMOUNT, DETREND_FRAMES);
    [cROIs, Cs] = mergeROIs3(Cs3, cROIs3, temporalCorrThr, spatialOverlapThr, fullSize, normalize_option, activity_option, SM_AMOUNT, DETREND_FRAMES);
   
    %% Clean up ROIs - Threshold above % of max, connected components above a minimum size
    options.thr_method = 'max';
    options.maxthr = max_pct;
    options.clos_op = strel('square',3);
    options.conn_comp = false;
    options.block_size = fullSize;
    options.min_pixel = min_sz;
    options.d1 = fullSize(1);
    options.d2 = fullSize(2);
    options.d3 = 1;
    cROIs = threshold_components3(cROIs,options);
    [cROIs,Cs] = component_split(cROIs,Cs,options.block_size,options.min_pixel,Inf);
    cROIs = sparse(cROIs);
    
    fprintf('Merging complete!\n');
end