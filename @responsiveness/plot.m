function [obj, varargout] = plot(obj,varargin)
%@responsiveness/plot Plot function for responsiveness object.
%   OBJ = plot(OBJ) 
%  The default plot without any optional input arguments is the
%  distribution of the responsiveness values either individually or for the
%  group object
%  The heterogeneity option is for pairs of cells within groups
%  To plot the absolute difference between all pairwise differences within
%  a group Ex. plot(obj,'Heterogeneity','IntraGroup','Pairs','Diff')
%
%  To plot the heterogeneity of the responsiveness values for each group
%  Ex. plot(obj,'Heterogeneity','IntraGroup','Bar')

Args = struct('showTitle',1, ...
    'Color','b','xlabel',0,'linkedZoom',0,'numBins',25, ...
    'NormalizeYaxis',0, 'Heterogeneity',0,'IntraGroup',0,...
    'InterGroup',0, 'Pairs',0,'STD',0,'Diff',0,'Mean',0,'Bar',0);  %yaxis probability rather than num of occurences

Args.flags = {'showTitle','xlabel','linkedZoom',...
    'Heterogeneity','IntraGroup',...
    'InterGroup','Pairs','STD','Diff','Mean','Bar'};

[Args,varargin2] = getOptArgs(varargin,Args,...
     'Color','numBins');

if Args.Heterogeneity  
    if Args.IntraGroup
        if Args.Pairs
            if Args.STD
                H = heterogeneity(obj,'IntraGroup','Pairs','STD');
            end
            if Args.Diff
                if Args.Mean
                    H = heterogeneity(obj,'IntraGroup','Pairs','Diff','Mean');
                else
                    H = heterogeneity(obj,'IntraGroup','Pairs','Diff');
                end
            end            
        elseif Args.Bar
            H = heterogeneity(obj,'IntraGroup','Bar'); 
        else
            if Args.STD
                H = heterogeneity(obj,'IntraGroup','STD');        
            end
            if Args.Diff
                H = heterogeneity(obj,'IntraGroup','Diff');
            end
        end 
        if Args.Bar
            resp_values=[];
            group_loc=0;
            rem=0;
            NumCells=[];
            for ii = 1:size(H,2)
                index = find(isfinite(H(:,ii))==1);
                if length(index) > 1  %%% Only include the responsiveness values for groups of cells
                    resp_values = [resp_values 0 0 flipdim(sort(H(index,ii)),1)'];   
                    group_loc = [group_loc group_loc(end)+2+(length(index)/2)+rem];
                    rem = length(index)/2;
                    NumCells = [NumCells length(index)];
                end
            end
            NumCells=sum(NumCells);
            group_loc(1)=[];
            ind = find(resp_values>=50);low_vals=resp_values;low_vals(ind)=0;
            bar(low_vals,1)
            set(gca,'FontSize',16)
            h = findobj(gca,'Type','patch');set(h,'FaceColor','r');hold on
            ind = find(resp_values<50);high_vals=resp_values;high_vals(ind)=0;
            bar(high_vals,1)   
            set(gca,'FontSize',16)
            set(gca,'XTick',group_loc)
            set(gca,'XTickLabel',[1:1:length(group_loc)])
        else
            binedges = 0:(100/Args.numBins):100; % Responsiveness Values range from 0 to 100.
            n = histcie(H,binedges);    
            bar(binedges,n,1,Args.Color)
            set(gca,'FontSize',16)
            NumCells=length(H);
        end
        if Args.Pairs
            if Args.STD
                ylabel('Number of Cell Pairs')
                xlabel('STD of Responsiveness Values')
            end
            if Args.Diff
                if Args.Mean
                    ylabel('Number of Groups')
                    xlabel('Mean Differences Between Group Responsiveness Values')
                else
                    ylabel('Number of Cell Pairs')
                    xlabel('Absolute Differences Between Responsiveness Values')
                end
            end
        elseif Args.Bar
            ylabel('RESPONSIVENESS')
            xlabel('Group Numbers')
        else            
            if Args.STD
                ylabel('Number of Groups')
                xlabel('STD of Responsiveness Values')
            end
            if Args.Diff
                ylabel('Number of Groups')
                xlabel('Absolute Differences Between Responsiveness Values')
            end
        end        
    elseif Args.InterGroup
        H = heterogeneity(obj,'InterGroup')
        binedges = 0:(100/Args.numBins):100; % Responsiveness Values range from 0 to 100.
        n = histcie(H,binedges);    
        bar(binedges,n,Args.Color) 
        set(gca,'FontSize',16)
        ylabel('Number of Groups')
        xlabel('STD of Responsiveness Values')
    end   
    
else 
    
    binedges = 0:(100/Args.numBins):100; % Responsiveness Values range from 0 to 100.
    n = histc(obj.data.responsiveness,binedges);
    
    if Args.NormalizeYaxis
        n = n/sum(n);
    end
    
    bar(binedges,n,1,Args.Color)
    set(gca,'FontSize',16)
    
    if Args.NormalizeYaxis
        ylabel('Probability')
        axis tight
        xlim([min(diff(binedges))*-1 min(diff(binedges))+100])
        xlabel('Responsiveness')
    else
        ylabel('Number of Occurences')  
        axis tight
        xlim([min(diff(binedges))*-1 min(diff(binedges))+100])
        xlabel('Responsiveness')
    end
    NumCells=get(obj,'Number');
    
end 

if size(get(obj,'SessionDirs'),2)==1
    SessionDirs = get(obj,'SessionDirs');
    string = [SessionDirs{1}];
else
    if Args.Pairs
        string = ['RESPONSIVENESS -- ' num2str(size(H,2)) ' Pairs'];
    else  
        string = ['RESPONSIVENESS -- ' num2str(NumCells) ' Cells'];
    end
end  

if Args.showTitle
    h = title(string);
    f = findobj(h,'type','text');
    set(f,'interpreter','none')
end


