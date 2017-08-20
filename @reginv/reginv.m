function obj = reginv(varargin)
%@reginv Constructor function for reginv class

Args = struct('RedoLevels',0,'SaveLevels',0,'Auto',0,...
    'NumFramesBack',5,'Lambda',[2 3 4 5],'NumPixels',40,...
    'PSTH',0,'Reps',0,'Test',0,'TestType','Convolve',...
    'FindRevCorrCenter',0,'FindHandMapCenter',0);

Args.flags = {'Auto','PSTH','Test','FindRevCorrCenter',...
        'FindHandMapCenter'};

[Args,modvarargin] = getOptArgs(varargin,Args, ...
    'subtract',{'RedoLevels','SaveLevels'}, ...
    'shortcuts',{'redo',{'RedoLevels',1}; 'save',{'SaveLevels',1}}, ...
    'remove',{'Auto'});

% variable specific to this class. Store in Args so they can be easily
% passed to createObject and createEmptyObject
Args.classname = 'reginv';
Args.matname = [Args.classname '.mat'];
Args.matvarname = 'ri';

numArgin = nargin;
if(numArgin==0)
    % create empty object
    obj = createEmptyObject(Args);
elseif( (numArgin==1) & isa(varargin{1},Args.classname))
    obj = varargin{1};
else
    % create object using arguments
    if(Args.Auto)
        % change to the proper directory
        % [pdir,cdir] = getDataDirs('session','relative','CDNow');
        % check for saved object
        if(ispresent(Args.matname,'file','CaseInsensitive') ...
                & (Args.RedoLevels==0))
            fprintf('Loading saved %s object...\n',Args.classname);
            l = load(Args.matname);
            obj = eval(['l.' Args.matvarname]);
        else
            % no saved object so we will try to create one
            % pass varargin in case createObject needs to instantiate
            % other objects that take optional input arguments
            obj = createObject(Args,modvarargin{:});
        end
        % change back to previous directory if necessary
        % if(isempty(cdir))
        %     cd(cdir)
        % end
    end
end

function obj = createObject(Args,varargin)
% Load the adjspikes and stiminfo files
cwd = pwd; sp = adjspikes('Auto'); cd ../..; st = stiminfo('Auto'); cd(cwd);
%%%% Get all the movie frame information %%%%%%%%%%%
x_off = (st.data.iniInfo.ScreenWidth - st.data.iniInfo.frame_cols)/2;
y_off = (st.data.iniInfo.ScreenHeight - st.data.iniInfo.frame_rows)/2;
frame_height = st.data.iniInfo.frame_rows;
frame_width = st.data.iniInfo.frame_cols;
moviename = st.data.iniInfo.stimulus_ext;
if strcmp(moviename,'.pgl') % Due to a change in the ini file
    moviename = st.data.iniInfo.stimulus_root;
end
screen_width = st.data.iniInfo.ScreenWidth;
screen_height = st.data.iniInfo.ScreenHeight;
movie_frame_numbers = repmat(st.data.iniInfo.start_frame:st.data.iniInfo.end_frame,1,st.data.catInfo.num_repetitions); 
if Args.FindRevCorrCenter
    %%%% Get the revcorr information %%%%%%%%%%%%%%
    rfwd=getEquivSession(varargin{:});
    if isempty(rfwd)
        obj = createEmptyObject(Args);
        return
    end
    cd(rfwd); cd ../..
    rfst = stiminfo('Auto');
    if ~isfield(rfst.data.iniInfo,'grid_center')
        rf_grid_center = [rfst.data.iniInfo.grid_x_center rfst.data.iniInfo.grid_y_center];
        grid_pixels = rfst.data.iniInfo.grid_x_size/rfst.data.iniInfo.m_seq_size(1);
    else
        rf_grid_center = rfst.data.iniInfo.grid_center';
        grid_pixels = rfst.data.iniInfo.grid_x_size/rfst.data.iniInfo.m_seq_size(1);
    end
    cd(rfwd);
    rev = revcorr('auto');
    rfs = rev.data.R;
    for i = 1:Args.NumFramesBack
        rf_image = imresize(rfs(:,:,i),grid_pixels);
        movie_image = zeros(screen_height,screen_width);movie_image(:)=NaN;    
        movie_image(rf_grid_center(2)-(size(rf_image,1)/2):rf_grid_center(2)+(size(rf_image,1)/2)-1,rf_grid_center(1)-(size(rf_image,2)/2):rf_grid_center(1)+(size(rf_image,2)/2)-1) = rf_image;
        rf_images(:,:,i) = movie_image;
    end
    [jnk,ndx]=max(abs(rf_images(:)-.5));
    [i1,i2,i3]=ind2sub(size(rf_images),ndx(1));
    center_point=[i1 i2];
    right_edge=400+(frame_width/2);
    left_edge=400-(frame_width/2);
    top_edge=300-(frame_height/2);
    bottom_edge=300+(frame_height/2);
    msize=size(rf_images,1);
    cp=center_point+Args.NumPixels;
    cn=center_point-Args.NumPixels;
    if (cp(1)>bottom_edge | cn(1)<top_edge | cp(2)>right_edge | cn(2)<left_edge)
        obj = createEmptyObject(Args);
        return
    end    
    clear rf_images rfs rfst rev
elseif Args.FindHandMapCenter
    center_point = [round(mean(st.data.iniInfo.ReceptiveField.CenterY)) round(mean(st.data.iniInfo.ReceptiveField.CenterX))];
else
    center_point = [300 400];
end

if Args.Test
    num_frames = st.data.catInfo.num_frames;
    movie_frame_numbers = movie_frame_numbers;
    Threshold = 30;
    [LL] = createreginvmatrix(Args.NumPixels);
    [R,SS] = testreginv(Args,Threshold,movie_frame_numbers,moviename,center_point,screen_height,screen_width,y_off,x_off,frame_height,frame_width,num_frames);
    for ii = 1:Args.NumFramesBack
        S = SS{ii};
        for l = 1:length(Args.Lambda)
            L = LL*(10^Args.Lambda(l)); % Laplacian Matrix
            M = [S;L]; 
            f = M\[R;zeros((Args.NumPixels^2),1)];
            data.reginv(:,:,ii,l) = reshape(f,Args.NumPixels,Args.NumPixels);
            data.lambda{l} = 10^Args.Lambda(l);
            data.R = R;
        end
        fprintf(['Frames Back ',num2str(ii) '\n'])
    end    
else
    %%%%%% Get the frame numberst with spikes %%%
    spikes = sp.data.adjSpiketrain;
    frames_vector = sp.data.adjFramePoints;
    if frames_vector(1) < 1 %% This is due to the conversion from 1 datapoint to ms
        frames_vector(1) = 0;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    spike_counts = histcie(spikes,frames_vector,'DropLast'); 
    if Args.PSTH
        num_frames = st.data.catInfo.num_frames;
        num_reps = st.data.catInfo.num_repetitions;
        sc = reshape(spike_counts,num_reps,num_frames);
        R = (mean(sc)/(st.data.catInfo.frame_duration/1000))';
        movie_frame_numbers = movie_frame_numbers(1:num_frames);    
    elseif Args.Reps
        num_frames = st.data.catInfo.num_frames;
        R = spike_counts(1:Args.Reps*num_frames);        
        movie_frame_numbers = movie_frame_numbers(1:Args.Reps*num_frames);
    else
        index = find(spike_counts>0);
        movie_frame_numbers = movie_frame_numbers(index); 
        R = spike_counts(index); 
    end
    clear st sp
    num_frames_back = (1:Args.NumFramesBack)-1;
    pix = Args.NumPixels/2;
    for ii = 1:Args.NumFramesBack
        movie_frame_numbers = movie_frame_numbers-num_frames_back(ii); 
        S = zeros(length(R),Args.NumPixels^2); %% Put into Memory the S Matrix
        tic
        for ff = 1:length(movie_frame_numbers)
            movie_frame = ReadPGLFrames(movie_frame_numbers(ff),moviename);
            screen_image = zeros(screen_height,screen_width);
            screen_image(y_off:frame_height+y_off-1,x_off:frame_width+x_off-1) = movie_frame;
            screen_image = screen_image(center_point(1)-(pix-1):center_point(1)+pix,center_point(2)-(pix-1):center_point(2)+pix);
            S(ff,:) = reshape(screen_image,1,Args.NumPixels^2);
        end   
        for l = 1:length(Args.Lambda)
            L = LL*(10^Args.Lambda(l)); % Laplacian Matrix
            M = [S;L]; 
            f = M\[R;zeros((Args.NumPixels^2),1)]; clear M L
            data.reginv(:,:,ii,l) = reshape(f,Args.NumPixels,Args.NumPixels);
            data.lambda{l} = 10^Args.Lambda(l);
            data.R = R;
        end  
        clear S
        toc
        fprintf(['Frames Back ',num2str(ii) '\n'])
    end
end

cd(cwd);
data.numSets = 1;
data.setNames{1} = cwd;
% create nptdata so we can inherit from it
n = nptdata(data.numSets,0,cwd);
d.data = data;
obj = class(d,Args.classname,n);
if(Args.SaveLevels)
    fprintf('Saving %s object...\n',Args.classname);
    eval([Args.matvarname ' = obj;']);
    % save object
    eval(['save ' Args.matname ' ' Args.matvarname]);
end	

function obj = createEmptyObject(Args)

% useful fields for most objects
data.numSets = 0;
data.setNames = '';

% these are object specific fields
data.reginv = [];
data.lambda = [];
data.R = [];

% create nptdata so we can inherit from it
n = nptdata(0,0);
d.data = data;
obj = class(d,Args.classname,n);