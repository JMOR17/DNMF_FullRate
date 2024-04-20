function [cROIs, Cs, coherence, skew, sz, patch_ID, A0] = DNMF_General5(V, options)
    % [cROIs, Cs, coherence, skew, sz] = DNMF_General3(V, options)
    
    %Extract options
    spatial_smooth = options.spatial_smooth;
    SNR_thr = options.SNR_thr;
    thr = options.thr;
    idealPatchSize = options.patchSize;
    stride = options.stride;
    overlapThr = options.overlapThr;
    sizeRange = options.sizeRange;
    filtSize = options.filtSize;
    DETREND_FRAMES = options.DETREND_FRAMES;
    shapeThr = options.shapeThr;
    med_opt = options.med_opt;
    height = size(V,1);
    width = size(V,2);

    CORE_OR_RANDOM = 1;
    %% Define ROI Pieces & Activity Traces
    
    % Split video into overlapping spatial patches
    
    indicesA = 1+(idealPatchSize(1):stride:height+stride-1)-idealPatchSize(1);
%     indicesA = 1:stride:(height-patchSize(1))+1;
    
    
    indicesB = 1+(idealPatchSize(2):stride:width+stride-1)-idealPatchSize(2);
%     indicesB = 1:stride:(width-patchSize(2))+1;
    
    ROIs = cell(1, length(indicesA)*length(indicesB));
    Cs = cell(length(indicesA)*length(indicesB), 1);
    COHERE = cell(length(indicesA)*length(indicesB), 1);
    SKEW = cell(length(indicesA)*length(indicesB), 1);
    SIZE = cell(length(indicesA)*length(indicesB), 1);
    PATCH_ID = cell(length(indicesA)*length(indicesB), 1);
    stamp = zeros(height, width);
    count = 1;
    % For each patch 
    for i_A = indicesA
        aa = i_A:min([(i_A+idealPatchSize(1)-1) height]);
        for i_B = indicesB
            bb = i_B:min([(i_B+idealPatchSize(2)-1) width]);
            fprintf('%d-%d\n',i_A,i_B);
            thisV0 = double(V(aa,bb,:));
            patchSize = [size(thisV0,1) size(thisV0,2)];
            if(spatial_smooth>0)
                thisV = gaussian_smooth_2d_slices(thisV0,spatial_smooth);
            else
                thisV = thisV0;
            end
            % Detrend this patch
            [dv,baseline] = j_detrend2b(thisV,DETREND_FRAMES,3,true,med_opt); 
    
            medThisV = nanmedian(dv.*zeroOut(dv>=0),3);
            dv(isnan(dv) | isinf(dv)) = 0;
            if(CORE_OR_RANDOM==1)
                % Find cores from original Y 
                A0 = findCores(dv, thr*medThisV, [1 1 filtSize],min(sizeRange));
                if(isempty(A0))
                    continue;
                end
                A0 = reshape(A0,[],size(A0,3));
                a0 = sparse(A0);
                roiAND = a0'*a0;
                roiOR = prod(patchSize) - (1-a0)'*(1-a0);
                roiJAC = roiAND./roiOR;
                [~,groups] = graphconncomp(roiJAC>overlapThr);
                temp = zeros(size(A0,1),max(groups));
                for i_group = 1:max(groups)
                    temp(:,i_group) = sum(A0(:,groups==i_group),2)>0;
                end
                A0 = temp;
                
%                 cThisV = reshape(thisV-baseline,[],size(thisV,3));
                cThisV = reshape(thisV,[],size(thisV,3));
                C0 = A0'*cThisV;
                SNR = quantile(C0,0.999,2)./mad(C0,[],2);
                valid = SNR>SNR_thr;
                A0 = A0(:,valid);
                C0 = C0(valid,:);
            else
                % or Initialize Randomly
                A0 = 50;
            end
            
            defoptions = CNMFSetParms;
            defoptions.d1 = size(thisV,1);
            defoptions.d2 = size(thisV,2);
            defoptions.block_size = patchSize;
            defoptions.min_pixel = sizeRange(1);
            defoptions.thr_method = options.thr_method;
            defoptions.maxthr = options.maxthr;
            defoptions.quantileThr = options.quantileThr;
            defoptions.final_C = options.final_C;
            defoptions.eta = options.eta;
            defoptions.beta = options.beta;
            if(isempty(A0))
                continue;
            end
            % Pass through CNMF for A, C estimates
            [A,C,B,F] = sparse_NMF_initialization7(thisV0,A0,defoptions);
            if(isempty(A))
                continue;
            end
            
            [coherence, skew, roiSize] = evaluateROIs(A, C, patchSize);
            
            % Calculate the correlation of the ROI with the video at the
            % frame of maximum activity, threshold based on that
            [~,i] = max(C,[],2);
            ppp = corrTrace(thisV0,reshape(full(A),[patchSize(1) patchSize(2) size(A,2)]),3);
            w = arrayfun(@(x) ppp(x,i(x)),1:size(i));
            valid = w>=shapeThr;
            
            A = A(:,valid);
            C = C(valid,:);
            coherence = coherence(valid);
            skew = skew(valid);
            roiSize = roiSize(valid);
            
            % Plotting
            temp = zeros(height, width, size(A,2));
            temp2 = zeros(height, width, size(A,2));           
            AA = reshape(bsxfun(@times,full(A),max(C,[],2)'), patchSize(1), patchSize(2), size(A,2));
            AA2 = reshape(full(A), patchSize(1), patchSize(2), size(A,2));
            
            temp(aa,bb,:) = AA;        
            temp2(aa,bb,:) = AA2;
            ROIs{count} = sparse(reshape(temp2,height*width,[]));
            Cs{count} = full(C);
            COHERE{count} = coherence(:);
            SKEW{count} = skew(:);
            SIZE{count} = roiSize(:);
            PATCH_ID{count} = count*ones(size(roiSize(:)));
            count = count+1;
            if(~isempty(temp))
                stamp = max(stamp, max(temp,[],3));
                clf;
                subplot(1,2,1);
                imagesc(stamp);
                axis square;
                subplot(2,2,2);
                imagesc(max(thisV,[],3));
                axis square;
                subplot(2,2,4);
                imagesc(reshape(max(A*C,[],2),[patchSize(1) patchSize(2)]));
                axis square;
                
                drawnow;
            end
        end
    end
    
    % Collect across patches
    fprintf('Collecting ROIs...');
    cROIs = cell2mat(ROIs);
    Cs = cell2mat(Cs);
    coherence = cell2mat(COHERE);
    skew = cell2mat(SKEW);
    sz = cell2mat(SIZE);
    patch_ID = cell2mat(PATCH_ID);
    fprintf('Done!\n');
end