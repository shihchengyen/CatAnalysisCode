function obj = plot(obj,varargin)
%@gsparsity/plot Plot function for gsparsity objects

Args = struct('Color',[0 0 1],'numBins',25,'Dist',0);    

Args.flags = {'Dist','SubPlot'};

[Args,modvarargin] = getOptArgs(varargin,Args);

if ~isempty(Args.NumericArguments)
    % plot one group at a time
    n = Args.NumericArguments{1};
    S = obj.data.groupsparsity(n);
    Values = obj.data.groupValues(:,n);
    Values = Values(find(Values>=0)); % Take out any NaN's
    name = obj.data.setNames{n};
else
    S = obj.data.groupsparsity;
    Values = obj.data.groupValues(:);
end

% Plot the Sparseness Values
if Args.Dist           
        [N] = flipud(sort(Values)); % Rank the stimulus counts in ascending order        
        plot(N,'Color',Args.Color,'LineStyle','-','Marker','*')
        title([name '  Group Sparsity: ' num2str(S)])
        xlabel('Stimulus Number')
        ylabel('Rate')      
else    
    binedges = 0:100/Args.numBins:100;
    [N] = histcie(S,binedges);
    bar(binedges,N,1)
    xlim([-((100/Args.numBins)/2) 100+((100/Args.numBins)/2)])
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor',Args.Color);
    set(gca,'FontSize',16)
    xlabel('Sparsity')
    ylabel('Number of Occurences')
    if exist('name')
        st = strfind(name,'all'); % Check for Groups or Pairs
        if ~isempty(st)
            title(['Sparsity of the Group ' name])
        else
            title(['Sparsity of the Pair ' name])
        end
    else
        st = strfind(obj.data.setNames{1},'all'); % Check for Groups or Pairs
        if ~isempty(st)
            title('Sparsity of the Groups')
        else
            title('Sparsity of Pairs within Groups')
        end
    end
end