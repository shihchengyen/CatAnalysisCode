function [obj,varargout] = plot(obj,varargin);
%latency/plot Plot function for LATENCY object

Args = struct('RefreshCorr',0,'AutoCorr',0,'FFT',0,'Color','b','Thresh',0.13, ...
    'RefreshesPerFrame',[85 3; 120 4; 150 6],'ReturnVars',{''},'SPTitleOff',0, ...
    'LabelsOff',0,'Replicate',0,'BarColor','b','LineColor',[1 0 0]);
Args.flags = {'RefreshCorr','AutoCorr','FFT','SPTitleOff','LabelsOff', ...
        'Replicate'};
Args = getOptArgs(varargin,Args);

if(isempty(Args.NumericArguments))
    % n not specified so plot all data in subplots
    numSets = get(obj,'Number');
    setVec = 1:numSets;
    if(Args.FFT)
        % can't just do nptFFTMag on obj.data.BinSums since they have
        % slightly different frequencies due to different number of points
        % for the different refresh rates
        % initialize population ratio variable
        pratio = zeros(numSets,1);
        framepratio = pratio;
        for idx = setVec
            % take the FFT
            [mx,f] = nptFFTMag(obj.data.BinSums{idx},1000);
            % find first f that is larger than the refresh rate
            refrate = obj.data.refreshRate(idx);
            ri = find(f>refrate);
            % take the mean of mx((ri-1):ri) and divide by mx(1) to get
            % ratio
            % pratio(idx) = mean(mx((ri(1)-1):ri(1)))/mx(1);
            pratio(idx) = mx(ri(1))/mx(1);
            % find number of refreshes per frame
            rpf = Args.RefreshesPerFrame;
            nframeref = rpf(find(rpf==round(refrate)),2);
            % find first f that is larger than the frame rate
            fi = find(f>refrate/nframeref);
            framepratio(idx) = mx(fi(1))/mx(1);
        end
        % find points above threshold
        pit = pratio > Args.Thresh;
        npit = ~pit;
        fit = framepratio > Args.Thresh;
        nfit = ~fit;
        % plot pratio
        plot(setVec(npit),pratio(npit),'.')
        hold on
        plot(setVec(pit),pratio(pit),'r.')
        plot(setVec(nfit),framepratio(nfit),'*')
        plot(setVec(fit),framepratio(fit),'r*')
        hold off
    else
        for i = 1:numSets
            binsums = obj.data.BinSums{i};
            sbinsize = obj.data.ActualSubBinSize(i);
            nsbins = length(binsums);
            nptSubplot(numSets,i);        
            % get xvals
            finalxval = (nsbins-1) * sbinsize;
            xvals = 0:sbinsize:finalxval;
            % plot distribution
            h = bar(xvals,binsums,'histc');
            set(h,'FaceColor',Args.BarColor,'EdgeColor',Args.BarColor);
            hold on
            % plot frequencies related to frequencies
            refresh = obj.data.refreshRate(i);
            refreshDuration = 1000/refresh;
            if(Args.Replicate)
                % replicate distribution to illustrate response latencies beyond the
                % frame limit
                h = bar(finalxval+xvals,binsums,'histc');
                set(h,'FaceColor','c');
                % get axis limits
                ax1 = axis;
                if(refresh>130)
                    % 150 Hz so there are 6 refreshes
                    line(repmat(refreshDuration * [1:12],2,1),repmat([ax1(3);ax1(4)],1,12), ...
                        'Color',Args.LineColor)
                elseif(refresh>90)
                    % 120 Hz so there are 4 refreshs
                    line(repmat(refreshDuration * [1:8],2,1),repmat([ax1(3);ax1(4)],1,8), ...
                        'Color',Args.LineColor)
                else
                    % 85 Hz so there are 3 refreshs
                    line(repmat(refreshDuration * [1:6],2,1),repmat([ax1(3);ax1(4)],1,6), ...
                        'Color',Args.LineColor)
                end
                xlim([0 finalxval+xvals(end)])
            else
                ax1 = axis;
                if(refresh>130)
                    % 150 Hz so there are 6 refreshes
                    line(repmat(refreshDuration * [1:6],2,1),repmat([ax1(3);ax1(4)],1,6), ...
                        'Color',Args.LineColor)
                elseif(refresh>90)
                    % 120 Hz so there are 4 refreshs
                    line(repmat(refreshDuration * [1:4],2,1),repmat([ax1(3);ax1(4)],1,4), ...
                        'Color',Args.LineColor)
                else
                    % 85 Hz so there are 3 refreshs
                    line(repmat(refreshDuration * [1:3],2,1),repmat([ax1(3);ax1(4)],1,3), ...
                        'Color',Args.LineColor)
                end
                xlim([0 finalxval])
            end
            hold off
            if(Args.SPTitleOff)
                title(num2str(i))
            else
                % compute FFT of the histograms to determine if the refresh
                % power exceed threshold
                [mx,f] = nptFFTMag(obj.data.BinSums{i},1000);
                % find first f that is larger than the refresh rate
                ri = find(f>obj.data.refreshRate(i));
                % take the mean of mx((ri-1):ri) and divide by mx(1) to get
                % ratio
                pratio = mx(ri(1))/mx(1);
                if(pratio>Args.Thresh)
                    rstring = '*';
                else
                    rstring = '';
                end
                title([getDataDirs('ShortName','DirString',obj.nptdata.SessionDirs{i}) rstring]);
            end
            ylim([0 max(binsums)])
            if(Args.LabelsOff)
                set(gca,'XTickLabel','','YTickLabel','')
            end
        end
        % select the bottom-left subplot
		nptSubplot(numSets,'BottomLeft');
        xlabel('Time (ms)')
        ylabel('Occurrences')
    end
else
    % plot for individual cell
    n = Args.NumericArguments{1};
    binsums = obj.data.BinSums{n};
    if(Args.AutoCorr)
        [c,lags] = xcorr(binsums,binsums);
        plot(lags,c)
    elseif(Args.RefreshCorr)
        % get delta functions at the refresh
        refresh = obj.data.refreshRate(n);
        refreshDuration = 1000/refresh;
        if(refresh>130)
            nrefresh = [1:6];
        elseif(refresh>90)
            nrefresh = [1:4];
        else
            nrefresh = [1:3];
        end
        refreshSignal = zeros(size(binsums));
        refreshSignal(round(nrefresh*refreshDuration)) = 1;
        [c,lags] = xcorr(binsums,refreshSignal);
        plot(lags,c)
    elseif(Args.FFT)
        [mx,f] = nptFFTMag(binsums,1000);
        plot(f,mx,'.-')
        hold on
        r = obj.data.refreshRate(n);
        % find first entry in f that is greater than r
        riplus = find(f>r);
        ri = riplus(1);
        % get the ratio of the mean of ri-1 and ri to DC
        % ratio = mean(mx((ri-1):ri))/mx(1);
        ratio = mx(ri)/mx(1);
        rpf = Args.RefreshesPerFrame;
        nframeref = rpf(find(rpf==round(r)),2);
        % find first f that is larger than the frame rate
        fiplus = find(f>r/nframeref);
        fi = fiplus(1);
        framepratio = mx(fi)/mx(1);
        ax1 = axis;
        line(repmat([r r/nframeref],2,1),repmat(ax1(3:4)',1,2),'Color','r')
        pind = [ri fi];
        plot(f(pind),mx(pind),'o')
        hold off
        zoom on
    else
        sbinsize = obj.data.ActualSubBinSize(n);
        nsbins = length(binsums);
        % get xvals
        finalxval = (nsbins-1) * sbinsize;
        xvals = 0:sbinsize:finalxval;
        % plot distribution
        h = bar(xvals,binsums,'histc');
        set(h,'FaceColor',Args.Color);
        hold on
        % replicate distribution to illustrate response latencies beyond the
        % frame limit
        h2 = bar(finalxval+xvals,binsums,'histc');
        set(h2,'FaceColor','c');
        % plot frequencies related to frequencies
        refresh = obj.data.refreshRate(n);
        refreshDuration = 1000/refresh;
        % get axis limits
        ax1 = axis;
        if(refresh>130)
            % 150 Hz so there are 6 refreshes
            line(repmat(refreshDuration * [1:12],2,1),repmat([ax1(3);ax1(4)],1,12), ...
                'Color','r')
        elseif(refresh>90)
            % 120 Hz so there are 4 refreshs
            line(repmat(refreshDuration * [1:8],2,1),repmat([ax1(3);ax1(4)],1,8), ...
                'Color','r')
        else
            % 85 Hz so there are 3 refreshs
            line(repmat(refreshDuration * [1:6],2,1),repmat([ax1(3);ax1(4)],1,6), ...
                'Color','r')
        end
        hold off
    end
    if(Args.FFT)
        title([getDataDirs('ShortName','DirString',obj.nptdata.SessionDirs{n}) ' Ratio: ' num2str(ratio)])
    else
        title(getDataDirs('ShortName','DirString',obj.nptdata.SessionDirs{n}));
    end
end

rvarl = length(Args.ReturnVars);
if(rvarl>0)
    % assign requested variables to varargout
    for rvi = 1:rvarl
        varargout{1}{rvi} = Args.ReturnVars{rvi};
        varargout{1}{rvi*2} = eval(Args.ReturnVars{rvi});
    end
end
