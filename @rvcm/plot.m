function [obj, varargout] = plot(obj,varargin)
%@rvcm/plot Plot function for rvcm object.

Args = struct('LabelsOff',0,'GroupPlots',1,'GroupPlotIndex',1,'Color','b', ...
		  'ReturnVars',{''}, 'ArgsOnly',0,'OverPlot',0,'ShowTitle',0, ...
          'Diff',0,'MeanStd',0);

Args.flags = {'LabelsOff','ArgsOnly','OverPlot','ShowTitle','Diff', ...
    'MeanStd'};

Args = getOptArgs(varargin,Args);

% if user select 'ArgsOnly', return only Args structure for an empty object
if Args.ArgsOnly
    varargout{1} = {'Args',Args};
    return;
end

if(isempty(Args.NumericArguments))
    % population plot
    num_cells=obj.data.numSets;
    colors = nptDefaultColors(1:num_cells);
    if(Args.MeanStd)
        d = obj.data.rvcmcn';
        maxreps = size(obj.data.rvcmcn,1);
        ha = errorbar(1:maxreps,nanmean(d),nanstd(d)./sqrt(sum(~isnan(d))));
        hb = get(ha,'Children');
        set(hb(1),'Color','r')
        Xdata = get(hb(2),'Xdata');
        temp = 4:3:length(Xdata);
        temp(3:3:end) = [];
        xleft = temp; 
        xright = temp+1;
        Xdata(xleft) = Xdata(xleft) + 1.5;
        Xdata(xright) = Xdata(xright) - 1.5;
        set(hb(2),'Xdata',Xdata);
        xlim([0 maxreps])
        line([0 maxreps],[0 0],'Color','g')
        xlabel('Repetition Number')
        ylabel('Change in Correlation Coefficent')
    else
        if(Args.Diff)
            for a = 1:num_cells
                plot(obj.data.rvcmcn(:,a),'.-','Color',colors(a,:));
                hold on
            end
        else
            for a = 1:num_cells
                plot(obj.data.rvcmc(:,a),'.-','Color',colors(a,:));
                [mxv,mxi]=max(obj.data.rvcmc(:,a));
                [mnv,mni]=min(obj.data.rvcmc(:,a));
                hold on
                plot(mxi,mxv,'^-','Color',colors(:,a),'LineWidth',2);
                plot(mni,mnv,'v-','Color',colors(:,a),'LineWidth',2);
            end
        end
        xlabel('Repetition Number')
        ylabel('Correlation Coefficent')
    end
else
    n = Args.NumericArguments{1};
    if(Args.Diff)
        plot(obj.data.rvcmcn(:,n));
        if(~Args.LabelsOff)
            ylabel('Change in Correlation Coefficent')
        end
    else
        plot(obj.data.rvcmc(:,n));
        [mxv,mxi]=max(obj.data.rvcmc(:,n));
        [mnv,mni]=min(obj.data.rvcmc(:,n));
        hold on
        plot(mxi,mxv,'^-','LineWidth',2);
        plot(mni,mnv,'v-','LineWidth',2);
        if(~Args.LabelsOff)
            ylabel('Correlation Coefficent')
        end
    end
    if(~Args.LabelsOff)
        xlabel('Repetition Number')
    end
    if Args.ShowTitle
        title(obj.data.setNames{n})
    end        
    hold off

end

rvarl = length(Args.ReturnVars);
if(rvarl>0)
     % assign requested variables to varargout
     for rvi = 1:rvarl
     	 rvidx = rvi * 2;
         varargout{1}{rvidx-1} = Args.ReturnVars{rvi};
         varargout{1}{rvidx} = eval(Args.ReturnVars{rvi});
     end
end

return
if(Args.Fast)
    raster = (obj.data.raster(:,(obj.data.rasterIndex(n)+1):obj.data.rasterIndex(n+1)));
    [rrows,rcols] = size(raster);
    % create the appropriate y offsets
    if(Args.Interleaved)
        ii = repmat(((1:rcols)+((Args.GroupPlotIndex-1)/(Args.GroupPlots+1))), ...
            rrows,1);
    else
        ii = repmat(1:rcols,rrows,1);
    end
    %%%%%% Start of Square Grating Hack Section %%%%%%%%%%%%
    % Load the stimInfo file for Square Grating Sessions and reorder the
    % original raster based on the orientations, temporal and spatial
    % frequencies and then plot a polar histogram for each spatial and
    % temporal frequency.
    sn=[obj.data.setNames{1}(9:10) obj.data.setNames{1}(26:27) '_stiminfo'];
    load(sn)
    if strcmpi(stimInfo.data.iniInfo.type,'sparse_noise')
        num_ori=stimInfo.data.catInfo.num_orientations;
        num_dir=stimInfo.data.catInfo.directions_per_orientation;
        num_sf=stimInfo.data.catInfo.num_spatial_frequencies;
        num_tf=stimInfo.data.catInfo.num_temporal_frequencies;
        num_reps=stimInfo.data.catInfo.num_repetitions;
        num_stim=num_ori*num_dir*num_sf*num_tf;
        seq=stimInfo.data.catInfo.seq_file(1:num_stim,:)+1; seq(:,1:2)=[]; seq(:,2)=[];
        nseq=[];
        for a=1:num_ori*num_dir
            o_ind=find(seq(:,1)==a);
            for b=1:num_tf
                t_ind=find(seq(:,2)==b);
                for c=1:num_sf
                    s_ind=find(seq(:,3)==c);
                    ind=intersect(o_ind,t_ind);
                    nseq=[nseq; intersect(ind,s_ind)];
                end
            end
        end
        oseq=nseq;
        for d=1:num_reps-1
            nseq=[nseq oseq+(num_stim*d)];
        end
        nseq=nseq';
        nseq=reshape(nseq,num_stim*num_reps,1);
        raster=raster(:,nseq);
        %% New order, Orientations, Temporal Frequency, Spatial Frequency

        % Create the polar histogram for each Spatial and Temporal Frequency
        % based on the orientation
        th=stimInfo.data.iniInfo.orientation_angles * (pi/180); %theta in radians for the polar plot
        sf=stimInfo.data.iniInfo.spatial_frequencies;
        tf=stimInfo.data.iniInfo.temporal_frequencies; tfi=1;
        bs=stimInfo.data.iniInfo.stimulus_duration_in_milliseconds;
        for e=1:num_tf*num_sf
            oind=((e*num_reps)-num_reps)+1;
            oind=oind:oind+(num_reps-1);
            for f=1:num_ori*num_dir
                data=raster(:,oind)';
                oind=oind+(num_tf*num_sf*num_reps);
                [y] = slidingHist(data,bs,bs,bs);
                N(f)=sum(y);
            end
            NN(e,:)=N;
            %             subplot(num_tf,num_sf,e)
            %             polar(th,N,'-')
            %             if e <= num_sf
            %                 title(num2str(sf(e)))
            %             end
            %             if e == ((tfi*num_sf)-num_sf)+1
            %                 ylabel(num2str(tf(tfi)))
            %                 tfi=tfi+1;
            %             end
        end
        [mv,mi]=max2(NN);
        NN=NN./mv;
        for e=1:num_tf*num_sf
            subplot(num_tf,num_sf,e)
            polar(0,1,'--k')
            hold on
            polar(th,NN(e,:),'-*')
            if e==mi
                polar(th,NN(e,:),'-*r')
            end
            if e <= num_sf
                if e==1
                    title(['SF: ',num2str(sf(e))])
                else
                    title(num2str(sf(e)))
                end
            end
            if e == ((tfi*num_sf)-num_sf)+1
                if e==1
                    ylabel(['TF: ',num2str(sf(e))])
                else
                    ylabel(num2str(tf(tfi)))
                end
                tfi=tfi+1;
            end
        end
        figure
    end
    %%%%%% End of Square Grating Hack Section %%%%%%%%%%%%
    if(isempty(Args.Highlight))
        plot(raster,ii,'LineStyle','none','Marker','.','Color',Args.Color);
    else
        hlight = Args.Highlight{Args.GroupPlotIndex};
        % figure out how many groups to highlight
        hgroups = length(hlight);
        % get the default marker oder
        dmarkers = nptDefaultMarkers(varargin{:});
        % save the first marker for the final plot
        marker1 = dmarkers(1);
        markerorder = dmarkers(2:end);
        % create logical array to keep track of which points have not been
        % plotted
        finalidx = logical(zeros(rrows,rcols));
        for hidx = 1:hgroups
            % get the indices for group hidx
            idx = hlight{hidx};
            % make sure idx is not larger than raster
            if(size(idx,1)>rrows)
                idx = idx(1:rrows,:);
            end
            % plot highlighted data
            plot(raster(idx),ii(idx), ...
                'LineStyle','none','Marker',markerorder(hidx),'Color',Args.Color);
            hold on
            % save the indices so we can plot any remaining un-highlighted
            % points at the end
            finalidx = finalidx | idx;
        end
        % get the indices for the rest of the data
        finalidx2 = ~finalidx;
        % plot the rest of the data
        plot(raster(finalidx2),ii(finalidx2), ...
            'LineStyle','none','Marker',marker1,'Color',Args.Color);
        hold off
    end
    ylim([-1 rcols+2])
    set(gca,'FontSize',12)
    ylabel('Repetitions','FontSize',12)
    xlabel('Time (ms)','FontSize',12)
    if(~Args.HideTitle)
        title(pwd)
    end
else
    % get relevant directory
    % dirnames = get(obj,'SessionDirs');
    % thisdir = dirnames{n};
    thisdir = obj.nptdata.SessionDirs{n};
    % save current directory
    cwd = pwd;
    % change to relevant directory
    cd(thisdir);

    if isempty(Args.chunkSize)
        sp=ispikes('auto');
        chunkSize = sp.data.chunkSize;
    else
        chunkSize = Args.chunkSize;
    end

    stimInfo = stiminfo('auto');
    time = [(n-1)*chunkSize ...
        n*chunkSize]*1000;
    rd = stimInfo.data.catInfo.repetition_duration;
    x=rem(time,rd);
    y=floor(time/rd);

    if Args.Stacked
        colors = 'brcgkmbrcgkmbrcgkm';
        ii = transpose([1:size(obj.data.raster',1)])*ones(1,size(obj.data.raster',2));
        plot(obj.data.raster',ii,'LineStyle','none','Marker','.','Color',Args.Color);
        group_num=dirnames{1}(37);
        for i = 1:obj.data.numSets
            hold on
            if dirnames{i}(37) ~= group_num
                line([0 stimInfo.data.catInfo.frame_duration*stimInfo.data.catInfo.num_frames],[obj.data.rasterIndex(i+1)+.5 obj.data.rasterIndex(i+1)+.5],'Color','r','LineWidth',2)
            else
                line([0 stimInfo.data.catInfo.frame_duration*stimInfo.data.catInfo.num_frames],[obj.data.rasterIndex(i+1)+.5 obj.data.rasterIndex(i+1)+.5],'Color','k','LineWidth',2)
            end
            group_num = dirnames{i}(37);
        end
        if(~Args.FMarkersOff)
            FrameMarkers = 0:stimInfo.data.catInfo.frame_duration:stimInfo.data.catInfo.frame_duration*stimInfo.data.catInfo.num_frames;
            plot(FrameMarkers,0,'^k','LineWidth',2)
            plot(FrameMarkers,obj.data.rasterIndex(end)+1,'vk','LineWidth',2)
        end
        ylim([0 obj.data.rasterIndex(end)+1])
    else
        %Create Raster
        num_repetitions = stimInfo.data.catInfo.num_repetitions;
        num_frames = stimInfo.data.catInfo.num_frames;
        frame_duration = stimInfo.data.catInfo.frame_duration;
        raster = (obj.data.raster(:,(obj.data.rasterIndex(n)+1):obj.data.rasterIndex(n+1)))';
        %Entire Raster (SLOW)
        if Args.showChunk==0   %if n==1 & Args.showChunk==0
            if size(obj.data.raster,2) > num_repetitions
                hold off %% this is for the displaydata gui
            end
            if Args.RefreshMarkers
                refreshMarks = 0:stimInfo.data.catInfo.sync_duration:frame_duration*num_frames;
                for i = 1:length(refreshMarks)
                    plot([refreshMarks(i) refreshMarks(i)],[0 num_repetitions+1],'--r','LineWidth',.5)
                    hold on
                end
            end
            if ~isempty(raster)
                if(Args.Interleaved)
                    [rrows,rcols] = size(raster);
                    ii = repmat(((1:rrows)+((Args.GroupPlotIndex-1)/Args.GroupPlots))', ...
                        1,rcols);
                else
                    ii = transpose([1:size(raster,1)])*ones(1,size(raster,2));
                end
                plot(raster,ii,'LineStyle','none','Marker','.','Color',Args.Color);
                if Args.Burst
                    burstraster = (obj.data.BurstRaster(:,obj.data.rasterIndex(n)+1:obj.data.rasterIndex(n)+num_repetitions))';
                    ii = transpose([1:size(burstraster,1)])*ones(1,size(burstraster,2));
                    hold on
                    plot(burstraster,ii,'LineStyle','none','Marker','.','Color','r');
                    hold off
                end
                set(gca,'YDir','reverse','FontSize',12)
            end
            %axis ij
            ylim([-1 num_repetitions+2])
            ylabel('Repetitions','FontSize',12)
            xlabel('Time','FontSize',12)
            if ~Args.FMarkersOff
                hold on
                FrameMarkers = 0:frame_duration:frame_duration*num_frames;
                plot(FrameMarkers,0,'vk','LineWidth',2)
                plot(FrameMarkers,num_repetitions+1,'^k','LineWidth',2)
            end
            xlim([0 frame_duration*num_frames+1])
            % instantiate responsiveness object
            if(~Args.HideTitle)
                resp = responsiveness('auto');
                title([pwd ' R = ' num2str(resp.data.responsiveness)])
            end
            hold off
        elseif Args.showChunk %Raster Chunk
            if ~isempty(raster.spikes)
                if x(2)>x(1)
                    ind = find((raster>x(1)) & (raster<x(2)));
                    [ii,jj] = ind2sub(size(raster),ind);
                    plot(raster(ind),ii,'LineStyle','none','Marker','.','Color',Args.Color);
                else
                    ind = find((raster>x(1)) & (raster<rd));
                    [ii,jj] = ind2sub(size(raster),ind);
                    plot(raster(ind),ii,'LineStyle','none','Marker','.','Color',Args.Color);

                    ind = find((raster.spikes>0) & (raster.spikes<x(2)));
                    [ii,jj] = ind2sub(size(raster),ind);
                    plot(raster(ind)+rd,ii,'LineStyle','none','Marker','.','Color',Args.Color);
                end
            end
            hold on
            if x(2)>x(1)
                axis([x(1) x(2) 0 num_repetitions+1]);
                patch([x(1) x(1) x(2) x(2)],[y(1)+.5  y(1)+1.5  y(1)+1.5  y(1)+.5],'b', 'FaceAlpha', .01);
            else
                axis([x(1) rd+x(2) 0 num_repetitions+1]);
                patch([x(1) x(1) rd rd],[y(1)+.5  y(1)+1.5  y(1)+1.5  y(1)+.5],'b', 'FaceAlpha', .01);
                patch([x(2)+rd x(2)+rd rd rd],[y(2)+.5  y(2)+1.5  y(2)+1.5  y(2)+.5],'b', 'FaceAlpha', .01);
            end
            hold off
            % axis ij
            ylabel('Repetitions')
            % xlabel('Time (msec)')
        end
        %hold off
    end

    % return to previous directory
    cd(cwd)
end
rvarl = length(Args.ReturnVars);
if(rvarl>0)
    % assign requested variables to varargout
    for rvi = 1:rvarl
        varargout{1}{rvi} = Args.ReturnVars{rvi};
        varargout{1}{rvi*2} = eval(Args.ReturnVars{rvi});
    end
end
