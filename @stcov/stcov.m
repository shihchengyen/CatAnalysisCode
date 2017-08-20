function obj = stcov(varargin)
%@stcov Constructor function for stcov class
%   OBJ = stcov('auto') attempts to create a stcov object by ...

Args = struct('RedoLevels',0,'SaveLevels',0,'Auto',0,'NumSpikes',0,'NumFramesBack',5,...
    'FindCenter',0,'NumPixels',10,'Interp',0,'NumInterp',2,'InterpType','nearest');                
Args.flags = {'Auto'};
[Args,modvarargin] = getOptArgs(varargin,Args, ...
    'subtract',{'RedoLevels','SaveLevels'}, ...
    'shortcuts',{'redo',{'RedoLevels',1}; 'save',{'SaveLevels',1}}, ...
    'remove',{'Auto'});

% variable specific to this class. Store in Args so they can be easily
% passed to createObject and createEmptyObject
Args.classname = 'stcov';
Args.matname = [Args.classname '.mat'];
Args.matvarname = 'stc';

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
cwd = pwd;
sp = adjspikes('Auto');
cd ../..
st = stiminfo('Auto');
cd(cwd);
if strcmp(st.data.iniInfo.type,'Movie')
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
        movie_image = zeros(screen_height,screen_width);    
        movie_image(rf_grid_center(2)-(size(rf_image,1)/2):rf_grid_center(2)+(size(rf_image,1)/2)-1,rf_grid_center(1)-(size(rf_image,2)/2):rf_grid_center(1)+(size(rf_image,2)/2)-1) = rf_image;
        rf_images(:,:,i) = movie_image;
    end
    if Args.FindCenter
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
    end
    
    %%%%%% Get the frames with spikes in them %%%%
    spikes = sp.data.adjSpiketrain;
    frames_vector = sp.data.adjFramePoints;
    if frames_vector(1) < 1 %% This is due to the conversion from 1 datapoint to ms
        frames_vector(1) = 0;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    spike_counts = histcie(spikes,frames_vector,'DropLast')';
    index = find(spike_counts>Args.NumSpikes);
    movie_frame_numbers = movie_frame_numbers(index);
    spike_counts = spike_counts(index);    
    num_frames_back = (1:Args.NumFramesBack)-1;
    if Args.Interp
        stim_size = size(framesexpand(zeros((Args.NumPixels*2)+1,(Args.NumPixels*2)+1),1/grid_pixels,Args.InterpType),2);
        stim_frames = zeros(stim_size,stim_size,length(movie_frame_numbers)); %% Put in Memory
    else
        stim_frames = zeros((Args.NumPixels*2)+1,(Args.NumPixels*2)+1,length(movie_frame_numbers)); %% Put in Memory
    end
    for f = 1:Args.NumFramesBack
        movie_frame_numbers = movie_frame_numbers-num_frames_back(f); 
        for ff = 1:length(movie_frame_numbers)
            movie_frame = ReadPGLFrames(movie_frame_numbers(ff),moviename);
            image_width = size(movie_frame,2);
            image_height = size(movie_frame,1);
            screen_image = zeros(screen_height,screen_width);
            screen_image(y_off:image_height+y_off-1,x_off:image_width+x_off-1) = movie_frame;
            screen_image = screen_image(center_point(1)-Args.NumPixels:center_point(1)+Args.NumPixels,center_point(2)-Args.NumPixels:center_point(2)+Args.NumPixels);
            if Args.Interp
                stim_frames(:,:,ff) = framesexpand(screen_image,1/grid_pixels,Args.InterpType);
            else
                stim_frames(:,:,ff) = screen_image
            end
        end       
        [evals,eframes] = imcov(stim_frames);
        data.evals{f} = evals;
        data.eframes{f} = eframes;
        data.stim_frames{f} = stim_frames;
%         m_stim_frames = zeros(stim_size,stim_size,size(stim_frames,3));
%         for ff = 1:size(stim_frames,3)
%             m_stim_frames(:,:,ff) = stim_frames(:,:,ff)*spike_counts(ff);
%         end
        data.sta{f} = mean(stim_frames,3);
        stim_frames=[];
        m_stim_frames=[];
    end
    
else % for the mseq        
    if Args.FindCenter
        eval(['cd ' cwd]);
        rf = revcorr('Auto',varargin{:});
        [jnk,ndx]=max(abs(rf.data.R(:)-.5));
        [i1,i2,i3]=ind2sub(size(rf.data.R),ndx(1));
        center_point = [i1 i2];
        msize = size(rf.data.R,1);
        if msize < Args.NumPixels*2
            Args.NumPixels = Args.NumPixels/2;
        end
        cp=center_point+Args.NumPixels;
        cn=center_point-Args.NumPixels;
        if (cp(1)>msize | cp(2)>msize | cn(1)<1 | cn(2)<1)
            obj = createEmptyObject(Args);
            return
        end   
        
    end
    % create object data
    % this is a valid object
    % these are fields that are useful for most objects
    num_frames_back = (1:Args.NumFramesBack)-1;
    for f = 1:Args.NumFramesBack
        [mseq_frames,spike_counts] = getmseqframes(st,sp,num_frames_back(f),Args);
        if Args.FindCenter
            mseq_frames = mseq_frames(center_point(1)-Args.NumPixels:center_point(1)+Args.NumPixels,center_point(2)-Args.NumPixels:center_point(2)+Args.NumPixels,:);
        end
        if Args.Interp
            mseq_frames = framesexpand(mseq_frames,Args.NumInterp,Args.InterpType);
        end
        [evals,eframes] = imcov(mseq_frames);
        data.evals{f} = evals;
        data.eframes{f} = eframes;
        data.stim_frames{f} = mseq_frames;
        %m_stim_frames = zeros(size(mseq_frames,1),size(mseq_frames,2),size(mseq_frames,3));
        %for ff = 1:size(mseq_frames,3)
        %    m_stim_frames(:,:,ff) = mseq_frames(:,:,ff)*spike_counts(ff);
        %end
        data.sta{f} = mean(mseq_frames,3);
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
data.evals=[];
data.eframes=[];
data.stim_frames=[];
data.sta=[];

% create nptdata so we can inherit from it
n = nptdata(0,0);
d.data = data;
obj = class(d,Args.classname,n);
