function [obj varargout] = variability(varargin)
%@variability/variability Constructor function for VARIABILITY object
%   OBJ = variability(VARARGIN)

Args = struct('RedoLevels',0,'SaveLevels',0, ...
    'Auto',0,'CellName','','SurrogateFF','framesg', ...
    'SurrogateTRE','framesgTRE','FrameBins',10, ...
    'LastGratingBugDate',[2;20;2004], ...
    'GratingRepetitionSet',2, ...
    'MatSearchString','*g*c*.mat', ...
    'NumFrames',1, ...
    'NumSGFiles',10, ...
    'NumSGPerFile',100, ...
    'ArgsOnly',0);
Args.flags = {'Auto', 'ArgsOnly'};

[Args,varargin2] = getOptArgs(varargin,Args, ...
    'shortcuts',{'redo',{'RedoLevels',1};'save',{'SaveLevels',1}}, ...
    'subtract',{'RedoLevels','SaveLevels'});

% if user select 'ArgsOnly', return only Args structure for an empty object
if Args.ArgsOnly
    varargout{1} = {'Args',Args};
    obj = createEmptyObject;
    return;
end

if nargin==0
    % create empty object
    obj = createEmptyObject;
elseif( (nargin==1) & (isa(varargin{1},'variability')) )
    obj = varargin{1};
else
    % create object using arguments
    if(Args.Auto)
        % check for saved object
        if(ispresent('variability.mat','file','CaseInsensitive') ...
                & (Args.RedoLevels==0))
            fprintf('Loading saved variability object...\n');
            % load saved object and exit
            l = load('variability.mat');
            obj = l.va;
        else
            % try to instantiate adjspikes and stiminfo objects
            as = adjspikes('auto',varargin2{:});
            % get session directory
            sdir = getDataDirs('session','relative');
            % get current directory
            cwd = pwd;
            % go to session directory
            cd(sdir);
            st = stiminfo('auto',varargin2{:});
            % return to previous directory
            cd(cwd);
            if(isempty(as) || isempty(st))
                % no saved object so try to create one
                % try to find name of .mat file
                a = nptDir(Args.MatSearchString);
                if(isempty(a))
                    % create empty object
                    obj = createEmptyObject;
                else
                    mat = load(a(1).name);
                    % hack to work with current data files
                    if(strcmp(Args.MatSearchString,'data.mat'))
                        % get cell name based on directory path
                        cwd = pwd;
                        % get cell name from last part of path
                        [p,cname] = nptFileParts(cwd);
                        % get group name from next part of path
                        [p,gname] = nptFileParts(p);
                        % get session name from next part of path
                        [p,sname] = nptFileParts(p);
                        % get site name from next part of path
                        [p,tname] = nptFileParts(p);
                        % get cat name from next part of path
                        [p,aname] = nptFileParts(p);
                        % reconstruct to get cell name
                        Args.CellName = [aname tname sname gname cname];
                    else
                        [fpath,fname] = nptFileParts(a(1).name);
                        Args.CellName = fname;
                    end
                    obj = createObject(mat.data,Args);
                end
            else
                Args.CellName = pwd;
                mat.cell_info.adjusted_spiketrain = ...
                    as.data.adjSpiketrain;
                mat.stimulus_info = st.data.iniInfo;
                mat.stimulus_info.adjusted_frames_vector = ...
                    as.data.adjFramePoints;
                obj = createObject(mat,Args);
            end
        end
    else
        % hopefully we were given a data structure
        % create object using data structure
        obj = createObject(varargin{1},Args);
    end
end

function obj = createEmptyObject

n = nptdata(0,0);
d.data.cellname = '';
d.data.scmean = [];
d.data.scstd = [];
d.data.cellid = [];
d.data.fano = [];
d.data.fanoSurrPercent = [];
d.data.fanoSurrZScores = [];
d.data.framebins = [];
d.data.entropy = [];
d.data.entropySurrPercent = [];
d.data.entropySurrZScores = [];
obj = class(d,'variability',n);

function obj = createObject(data,Args)

% initialize flag
isGrating = 0;

index = 1;

d.data.cellname = {Args.CellName};
framebins = Args.FrameBins;

% compute spike counts in frames
datasc = histcie(data.cell_info.adjusted_spiketrain, ...
    data.stimulus_info.adjusted_frames_vector,'DropLast');

% check if this is a grating session
if(isfield(data.stimulus_info,'type') ...
        && isfield(data.stimulus_info,'obj_type') ...
        && strcmp(data.stimulus_info.type,'sparse_noise') ...
        && strcmp(data.stimulus_info.obj_type,'Grating'))
    % set flag so we don't have to check again
    isGrating = 1;
    % select stimulus to use and exclude ISI period
    % need a version number of some sort to distinguish between old
    % grating sessions that have bugs and future grating sessions
    % convert date to number since we only have date and not the
    % version of Presenter
    lastBugDate = getDateNumber(Args.LastGratingBugDate);
    thisDate = getDateNumber(data.stimulus_info.Date);
    if(thisDate<=lastBugDate)
        % compute number of grating frames presented
        % get number of repetitions
        nreps = data.stimulus_info.num_blocks;
        % reshape into nreps x stimulid.data.fano = (d.data.scstd.^2) ./ d.data.scmean;
        stimorder = reshape(data.stimulus_info.raster_order,nreps,[]);
        % get total number of stimuli
        nstimuli = nreps*length(data.stimulus_info.orientation_angles);
        % reshape sc into (num_frames+isi frame) x (numreps*stimuli)
        % need to add 0 at the end since in this version of Presenter
        % the last repetition is missing the ISI. adjusted_frames_vector
        % also contains an extra 0 at the front to account for the fact
        % that the first frame is missing on the first presentation
        screps = reshape(veccat(datasc,0),[],nstimuli);
        % index into screps to get the orientation we are interested in
        scori = screps(:,stimorder(:,Args.GratingRepetitionSet));
        % check if the first number in the column of stimorder that
        % we are interested is 1
        if(stimorder(1,Args.GratingRepetitionSet)==1)
            % first frame of first rep is missing so toss the data for
            % the first frame. Also toss the last frame since that is
            % the ISI period
            % we are going to take mean and std across columns so we
            % will just keep sc1 in num_frames x numreps
            datasc1 = scori(2:data.stimulus_info.num_frames,:);
            nframes = data.stimulus_info.num_frames - 1;
        else
            datasc1 = scori(1:data.stimulus_info.num_frames,:);
            nframes = data.stimulus_info.num_frames;
        end
    end
else
    % should be movie session
    % compute number of frames presented
    nframes = data.stimulus_info.end_frame ...
        - data.stimulus_info.start_frame + 1;
    % reshape into repetitions and then transpose so
    % we can take the mean and std easily
    datasc1 = reshape(datasc,nframes,[]);
    
end

repetitions = size(datasc1, 2);
aframes = Args.NumFrames;
numsgperfile = Args.NumSGPerFile;
nsgfiles = Args.NumSGFiles;
concatScmean = zeros(nframes-Args.NumFrames+1 ,numsgperfile*nsgfiles);
concatScstd = concatScmean;
jvec = 1:numsgperfile;
framebins = Args.FrameBins;
d.data.framebins = framebins;

if (aframes == 1)
    % keep track of cell id which will make it easier to identify cells
    % when variability objects are added together
    d.data.cellid = ones(nframes,1);
    d.data.scmean = mean(datasc1,2);
    d.data.scstd = std(datasc1,0,2);
    warning off MATLAB:divideByZero
    d.data.fano = (d.data.scstd.^2) ./ d.data.scmean;
    warning on MATLAB:divideByZero

    % check for combined FF file
    combinedFFname = [Args.SurrogateFF 'FF.mat'];
    combinedFF = nptDir(combinedFFname,'CaseInsensitive');
    if(~isempty(combinedFF))
        % load combined FF file
        mat = load(combinedFF.name);
        fn = fieldnames(mat);
        SData = getfield(mat,fn{1});
        concatScmean = SData.scmean;
        concatScstd = SData.scstd;
    else
        % no combined file found so load each individually
        for i = 1:nsgfiles
            idx = (((i-1)*numsgperfile+1):(i*numsgperfile));
            fname = [Args.SurrogateFF num2str(i) 'FF.mat'];
            if(ispresent(fname,'file'))
                load(fname);
                concatScmean(:,idx) = sg.scmean;
                concatScstd(:,idx)  = sg.scstd;
            end
        end
        % save combined file so we won't have to do this again
        SData.scmean = concatScmean;
        SData.scstd = concatScstd;
        save(combinedFFname,'SData');
    end

    % compute Fano for surrogates, dimensions are (numframes x surrsets)
    warning off MATLAB:divideByZero
    surrFano = (concatScstd.^2) ./ concatScmean;
    warning on MATLAB:divideByZero
    % get number of sets of surrogates
    surrSets = size(surrFano,2);
    % replicate data Fano
    dataFano = repmat(d.data.fano,1,surrSets);
    % compute percent of surrogates with higher Fano than data for each
    % frame
    d.data.fanoSurrPercent = sum(dataFano<surrFano,2) / surrSets;
    % get mean and std of surrogate data
    surrMean = mean(surrFano,2);
    surrStd = std(surrFano,0,2);
    warning off MATLAB:divideByZero
    d.data.fanoSurrZScores = (d.data.fano - surrMean) ./ surrStd;
    warning on MATLAB:divideByZero

    % compute entropy
    % check if this is a grating session
    if(isGrating)
        if(thisDate<=lastBugDate)
            % get bin limits for histcie
            afv = vecc(data.stimulus_info.adjusted_frames_vector);
            % get length of afv
            afvl = length(afv);
            dafv = diff(afv);
            binsizes = dafv/framebins;
            mat1 = tril(ones(framebins));
            mat2 = [afv(1:(afvl-1))'; repmat(binsizes',(framebins-1),1)];
            blimits = mat1 * mat2;
            binlimits = [reshape(blimits,[],1); afv(afvl)];
            sc = histcie(data.cell_info.adjusted_spiketrain, ...
                binlimits,'DropLast');
            % reshape into framebins x num_frames x nstimuli
            screps = reshape([sc;zeros(framebins,1)],framebins,[],nstimuli);
            % pick out the relevant trials from dimension 3
            scori = screps(:,:,stimorder(:,Args.GratingRepetitionSet));
            if(stimorder(1,Args.GratingRepetitionSet)==1)
                % first frame of first rep is missing so toss the data for the
                % first frame. Also toss the last frame since that is the ISI
                % period
                sc1a = scori(:,2:data.stimulus_info.num_frames,:);
            else
                sc1a = scori(:,1:data.stimulus_info.num_frames,:);
            end
            % reshape into framebins x (nframes * nreps)
            datasc1 = reshape(sc1a,framebins,[]);
        end
    else
        % get bin limits for histcie
        afv = vecc(data.stimulus_info.adjusted_frames_vector);
        % get length of afv
        afvl = length(afv);
        dafv = diff(afv);
        binsizes = dafv/framebins;
        mat1 = tril(ones(framebins));
        mat2 = [afv(1:(afvl-1))'; repmat(binsizes',(framebins-1),1)];
        blimits = mat1 * mat2;
        binlimits = [reshape(blimits,[],1); afv(afvl)];
        sc = histcie(data.cell_info.adjusted_spiketrain, ...
            binlimits,'DropLast');
        % reshape into framebins x (nframes * nreps)
        sc1 = reshape(sc,framebins,[]);
    end
    dentropy = getTREntropy(sc1,nframes);
    d.data.entropy = dentropy';
    trefname = [Args.SurrogateTRE '.mat'];
    if(ispresent(trefname,'file'))
        % load the data from the surrogates and compute percentage of surrogates
        % that have higher entropy
        stre = load(trefname);

        % surrogate entropies in matrix with dimensions surrsets x nframes
        surrSets = size(stre.entropy,1);
        % replicate data entropy
        dataentropy = repmat(dentropy,surrSets,1);
        %      dataentropy = repmat(d.data.entropy', surrSets, 1);
        % compute percent of surrogates with higher entropy than data for each
        % frame
        d.data.entropySurrPercent = (sum(dataentropy<stre.entropy) / surrSets)';
        % compute mean and std of surrogate data
        surrMean = mean(stre.entropy);
        surrStd = std(stre.entropy);
        warning off MATLAB:divideByZero
        d.data.entropySurrZScores = ((dentropy - surrMean) ./ surrStd)';
        warning on MATLAB:divideByZero
    else
        d.data.entropySurrPercent = zeros(nframes,1);
        d.data.entropySurrZScores = d.data.entropySurrPercent;
    end

else % (aframes > 1)
    %##################################################################
    surrFiles = Args.NumSGFiles;
    surrSetsPerFile = Args.NumSGPerFile;
    
    %##################################################################
    sbinsperframe = framebins/aframes;
    nfullframes = nframes - aframes + 1;
    ntotalframes = nfullframes * repetitions;

    d.data.cellid = ones(nframes-aframes+1,1);

    % 1. For DATA

    % 1.1 For dataFano calculation
    convVec = ones(1, Args.NumFrames);
    convDataSC = ((convmtx(convVec, nframes-Args.NumFrames+1)) * datasc1)';
    convDataScmean = mean(convDataSC)';
    convDataScstd = std(convDataSC)';
    d.data.scmean = convDataScmean;
    d.data.scstd = convDataScstd;
    d.data.fano = ((convDataScstd.^2) ./ convDataScmean );

    % 1.2 For dataEntropy
    afv = vecc(data.stimulus_info.adjusted_frames_vector);
    % get length of afv
    afvl = length(afv);
    if(~isGrating)

        % 1.2.1 For dataEntropy when sbinsperframe == 1
        if (rem(sbinsperframe,1)==0)
            dafv = diff(afv);
            binsizes = dafv/sbinsperframe;
            mat1 = tril(ones(sbinsperframe));
            mat2 = [afv(1:(afvl-1))'; repmat(binsizes',(sbinsperframe-1),1)];
            blimits = mat1 * mat2;
            binlimits = [blimits(:); afv(afvl)];
            % binlimits = [reshape(blimits, [], 1); afv(afv1)];
            sc = histcie(data.cell_info.adjusted_spiketrain, ...
                binlimits,'DropLast');
            % reshape into framebins x (nframes * nreps)
            i1 = reshape(0:(nframes*repetitions-1),nframes,[]);
            i2 = i1(1:nfullframes,:) * sbinsperframe + 1;
            i3 = reshape(i2,1,[]);
            i4 = tril(ones(framebins));
            i5 = [i3; ones(framebins-1,ntotalframes)];
            i6 = i4 * i5;
            sc1 = sc(i6);

            dentropy = getTREntropy(sc1, nframes-aframes+1);
            d.data.entropy = dentropy';

            % 1.2.2 For dataEntropy when sbinsperframe ~= 1
        else %(rem(sbinsperframe,1)~=0)
            frameNumber = zeros(aframes, 1);
            % Drop last point and reshape afv to nframes x repetitions
            afv2 = reshape(afv(1:afvl-1), nframes, repetitions);
            % Add first row and the last point to the last row
            afv3 = [afv2; afv2(1, 2:repetitions) afv(afvl)];
            % Loop tp get offsets
            for i = 1: aframes
                % Start from ith frame, and pick up the point every aframes
                afvtmp = afv3(i:aframes:nframes+1, :);
                % Record the number of new frames in each repettions
                frameNumber(i, :) = size(afvtmp, 1) - 1;
                % Reshape the new matrix into vector array
                afvtmp2 = [reshape(afvtmp, [], 1); afv(afvl)];
                % Take the difference
                dafvtmp2 = diff(afvtmp2);
                % Get binsizes
                binsizes = dafvtmp2/framebins;
                mat1 = tril(ones(framebins));
                mat2 = [afvtmp2(1:(length(afvtmp2)-1))'; repmat(binsizes', (framebins-1), 1)];
                blimits = mat1 * mat2;
                binlimits = [reshape(blimits, [], 1); afvtmp2(length(afvtmp2))];
                sc{i} = histcie(data.cell_info.adjusted_spiketrain, ...
                    binlimits,'DropLast');
                % Reshape histogram sc into repetitions x ...
                mat3 = reshape(cell2mat(sc(1, i)), repetitions, []);
                % Drop last column
                mat4 = mat3(:, 1:size(mat3, 2)-framebins);
                % Take transpose
                mat5 = mat4';
                % Reshape into framebins x ...
                sc{i} = reshape(mat5, framebins, []);
                afvtmp = [];
            end

            frameSum = sum(frameNumber(1:aframes));
            totalFrameNumber = frameSum*repetitions;
            vecA = 1:totalFrameNumber;
            vecB = reshape(vecA, [], repetitions);
            vecC = vecB';
            sc1 = zeros(framebins, totalFrameNumber);

            for idx = 1:aframes
                indaxCount = reshape(vecC(:, idx:aframes:frameSum)', 1, []);
                sc1(:, indaxCount) = cell2mat(sc(1, idx));
            end

            dentropy = getTREntropy(sc1, nframes-aframes+1);
            d.data.entropy = dentropy';
        end
    end
    clear sc;

    % 2. For SURROGATES

    % load refractory data
    mat = load('refractory.mat');
    rf = struct(mat.rf);

    % allocate memory
    sptrain = cell(1,surrFiles);
    convSurrScmean = zeros(nfullframes, surrFiles*surrSetsPerFile);
    convSurrScstd = convSurrScmean;
    surrEntropy = zeros(surrFiles*surrSetsPerFile, nframes-aframes+1);

    % get frame boundaries
    rtlength = length(rf.data.rtEdges{1});
    % take out last value in rtEdges so there will be an even number of points
    rt1 = reshape(rf.data.rtEdges{1}(1:(rtlength-1)),rf.data.qtframebins,[]);
    % get the last value which will be used a couple of times
    rtlast = rf.data.rtEdges{1}(rtlength);

    % Preparation for 2.1 For surrFano
    % grab the first row which should be the frame limits and add the last
    % point in rtEdges
    fbins = [rt1(1,:)'; rf.data.rtEdges{1}(rtlength)];

    % Preparation for 2.2.1 For surrEnttropy when rem(sbinsperframe,1)==0
    if (rem(sbinsperframe,1)==0)
        surrTotal = surrFiles * surrSetsPerFile;
        surrEntropy = zeros(numsgperfile*nsgfiles, nframes-aframes+1);
        % grab the first row which should be the frame limits
        frameLimits = [rt1(1,:) rtlast];
        % get difference between frame limits
        dfL = diff(frameLimits);
        % get bin sizes per frame
        binsize = dfL/sbinsperframe;
        % create matrix of bin limits
        mat1 = tril(ones(sbinsperframe));
        mat2 = [frameLimits(1:nframes); repmat(binsize,sbinsperframe-1,1)];
        binLimits = mat1 * mat2;
        % append time of last frame
        histEdges = [reshape(binLimits,[],1); rtlast];
        % Preparation for 2.2.1 For surrEnttropy when sbinsperframe ~= 1
    else % (rem(sbinsperframe,1)~=0)
        frameNumber = zeros(aframes, 1);
        % grab the first row which should be the frame limits
        frameLimits = [rt1(1,:) rtlast];
        % End of preparation
    end

    % Loop for calculate surrFano, surrEntropy when sbinsperframe==1 and
    % surrEntropy when sbinsperframe~=1
    % read surrogates
    % loop over surrogates
    
    for i = 1:surrFiles
        sptrain = readSurrogateBin([Args.SurrogateFF num2str(i) '.bin']);
        % loop over sets in each file
        for j = 1: surrSetsPerFile
            % get spike counts for each frame of sptrain, one column for each rep
            sc_for_fano = histcie(cell2array(sptrain{j}),fbins,'DropLast');
            convSurrSC = ((convmtx(convVec, nframes-aframes+1)) * sc_for_fano);
            % reshape into repetitions and then transpose so we can take the mean
            % and std easily
            convSurrSC1 = reshape(convSurrSC, nframes-aframes+1, [])';
            % store mean and std in columns
            concatSurrScmean(:,100*(i-1)+j) = mean(convSurrSC1)';
            concatSurrScstd(:,100*(i-1)+j) = std(convSurrSC1)';
            
            if (rem(sbinsperframe,1)==0)
                % binSC dimensions should be (numframes * framebins) x repetitions
                sc = histcie(cell2array(sptrain{j}),histEdges,'DropLast');
                % reshape into framebins x (nframes * nreps)
                i1 = reshape(0:(nframes*repetitions-1),nframes,[]);
                i2 = i1(1:nfullframes,:) * sbinsperframe + 1;
                i3 = reshape(i2,1,[]);
                i4 = tril(ones(framebins));
                i5 = [i3; ones(framebins-1,ntotalframes)];
                i6 = i4 * i5;
                sc1 = sc(i6);

                surrEntropy1 = getTREntropy(sc1,nframes-aframes+1);
                surrEntropy(100*(i-1)+j,:) = surrEntropy1;
            else %(rem(sbinsperframe,1)~=0)
                % Loop over aframes to get offset
                for k = 1:aframes
                    % Start from kth frame, and pick up the point every aframes
                    frameLimitsTmp = frameLimits(k:aframes:length(frameLimits));
                    % Get length of the new matrix frameLimitsTmp
                    frameLimitsTmpl = length(frameLimitsTmp);
                    % Take the difference
                    dframeLimitsTmp = diff(frameLimitsTmp);
                    % Record the number of new frames in each repettions
                    frameNumber(k, :) = size(frameLimitsTmp, 2) - 1;
                    % Get the binsizes
                    binsize = dframeLimitsTmp/framebins;
                    mat1 = tril(ones(framebins));
                    mat2 = [frameLimitsTmp(1:frameLimitsTmpl-1); repmat(binsize, (framebins-1), 1)];
                    blimits = mat1 * mat2;
                    binlimits = [reshape(blimits, [], 1); frameLimitsTmp(frameLimitsTmpl)];
                    sc{k} = histcie(cell2array(sptrain{j}),binlimits,'DropLast');
                    mat3 = reshape(cell2mat(sc(1,k)), 10, []);
                    sc{k} = mat3;
                    frameLimitsTmp = [];
                end
                frameSum = sum(frameNumber(1:aframes));
                totalFrameNumber = frameSum*repetitions;
                vecA = 1:totalFrameNumber;
                vecB = reshape(vecA, [], repetitions);
                vecC = vecB';
                sc1 = zeros(framebins, totalFrameNumber);

                for idx = 1:aframes
                    indaxCount = reshape(vecC(:, idx:aframes:frameSum)', 1, []);
                    sc1(:, indaxCount) = cell2mat(sc(1, idx));
                end

                surrEntropy(100*(i-1)+j,:) = getTREntropy(sc1, nframes-aframes+1);
            end
        end
    end

    warning off MATLAB:divideByZero
    surrFano = ((concatSurrScstd.^2) ./ concatSurrScmean);
    warning on MATLAB:divideByZero
    surrSets = size(surrFano, 2);
    dataFano = repmat(d.data.fano, 1, surrSets);
    d.data.fanoSurrPercent = sum(dataFano<surrFano, 2) / surrSets;
    surrMean = mean(surrFano,2);
    surrStd = std(surrFano,0,2);
    warning off MATLAB:divideByZero
    d.data.fanoSurrZScores = (d.data.fano - surrMean) ./ surrStd;
    warning on MATLAB:divideByZero

    % surrogate entropies in matrix with dimensions surrsets x nframes
    surrSets = size(surrEntropy,1);
    % replicate data entropy
    dataentropy = repmat(dentropy,surrSets,1);
    % compute percent of surrogates with higher entropy than data for each
    % frame
    d.data.entropySurrPercent = (sum(dataentropy<surrEntropy) / surrSets)';
    % compute mean and std of surrogate data
    surrMean = mean(surrEntropy);
    surrStd = std(surrEntropy);
    warning off MATLAB:divideByZero
    d.data.entropySurrZScores = ((dentropy - surrMean) ./ surrStd)';
    warning on MATLAB:divideByZero

end

% create objects
n = nptdata(1,0,'SessionDirs',{pwd});
obj = class(d,'variability',n);
if(Args.SaveLevels>0)
    fprintf('Saving variability object...\n');
    va = obj;
    save variability va
end

function n = getDateNumber(datevector)

% convert datevector to strings
d = num2str(datevector);
n = datenum([d(1,:) '/' d(2,:) '/' d(3,:)]);
