function [R,H,Surrogate_Inv] = FitSurrogateFR(files,st,Args)
%% This function is called in the firingrate plot function to calculate the
%% fit of the surrogate spiketrains to various functions

surrFiles = 1;%length(files);
surrSetsPerFile = 100;
histEdges = [0 st.data.catInfo.repetition_duration];
index = 1;
Surrogate_Inv=[];
% loop over surrogates
for i = 1:surrFiles
    % load surrogates
    sptrain = readSurrogateBin(['framesg' num2str(i) '.bin']);
    % loop over sets in each file
    for j = 1:surrSetsPerFile
        Raster = cell2array(sptrain{j});
        FR=[];
        if Args.Overlap % If a specified overlap is desired, the output will be a single vector
            
            m = round(Args.Binsize/Args.Overlap);
            shape = ones(1,m);
            edges = (histEdges(1):Args.Overlap:histEdges(end))';
            for r = 1:size(Raster,2)
                spikeCounts = histc(Raster(:,r),edges);
                spikeCounts = conv(spikeCounts,shape);
                spikeCounts(1:m-1)=[];spikeCounts(end-m:end)=[];
                FR = [FR;spikeCounts];
            end

        elseif Args.Count % If a specified SpikeRate is desired, the output will be a single vector

            mean_sc = length(find(Raster>=0))/((histEdges(end)*size(Raster,2))/1000);
            Args.Binsize = (Args.Count/mean_sc)*1000; % Get the binsize for a desired rate in spikes/second
            edges = (histEdges(1):Args.Binsize:histEdges(end))';
            for r = 1:size(Raster,2)
                spikeCounts = histcie(Raster(:,r),edges,'DropLast');
                FR = [FR;spikeCounts];
            end

        elseif Args.InstFR

            FR = 1./(diff(Raster)/1000);
            data=[]; %%% Reshape the data and take out the NaN's
            for r = 1:size(FR,2)
                non_zeros = find(FR(:,r)>=0);
                data = [data;FR(non_zeros,r)];
            end
            FR = data;
            %% Check for Inf for some reason %%
            ind = max(data);
            if ind == Inf
                ind = find(FR==ind);
                FR(ind) = -1;
                non_zeros = find(FR>=0);
                FR = FR(non_zeros);
                Surrogate_Inv = [Surrogate_Inv index];
            end
            
        else
            
            edges = (histEdges(1):Args.Binsize:histEdges(end))';
            for r = 1:size(Raster,2)
                spikeCounts = histcie(Raster(:,r),edges,'DropLast');
                FR = [FR;spikeCounts];
            end
            
        end

        bins = (min(FR(:)):(max(FR(:))-min(FR(:)))/Args.numBins:max(FR(:)))';
        N = histcie(FR,bins,'DropLast'); bins = bins(2:end);
        non_zeros = find(N>0);
        N = log10(N(non_zeros));
        bins = bins(non_zeros);
        H{index}(:,1) = bins;
        H{index}(:,2) = N;        
        %% Fit types to be used
        FitName1 = 'poly1';
        FitName2 = fittype('a + b*log(x)');
        %% Fit the data
        [fitobj,goodness,output] = fit(bins,N,FitName1); % Exp
        R(index,1) = goodness.adjrsquare;
        [fitobj,goodness,output] = fit(bins,N,FitName2); % Power
        R(index,2) = goodness.adjrsquare;
        index = index + 1;
        
    end
end
