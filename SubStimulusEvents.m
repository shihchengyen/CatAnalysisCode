% function all_ss_events = substimulusevents(binsize,fig,varargin);

binsize = 3;
fig = 1;

load jonsdisplayevents

bin = [3 5];
rez = find(bin == binsize);

ncells = length(displayevents);

% loop over cells
for e = 1:ncells
    
    raster = displayevents(e).raster;
    
    % find sub-stimulus events for cell e
    ss_event_ind = displayevents(e).ss_event_ind{rez};
            
    if isempty(ss_event_ind);
        % if there were no sub-stimulus events
        all_ss_events(e).ss_iei_mat = [];
    else
		% find limits of sub-stimulus events for cell e
		pos_lim = displayevents(e).first_pos_limit{rez};
		pos_lim1 = pos_lim(ss_event_ind);
		pos_lim2 = pos_lim(ss_event_ind+1);
		
		neg_lim = displayevents(e).first_neg_limit{rez};
		neg_lim1 = neg_lim(ss_event_ind);
		neg_lim2 = neg_lim(ss_event_ind+1);
		
		% get number of sub-stimulus events
		nss = length(pos_lim1);
		ss_events = zeros(1,nss);
		
		% loop over sub-stimulus events in cell e
		for x = 1:nss
			% get number of repetitions - can't move this outside for 
			% loop since some cells have different number of reps
			reps = length(raster);
			ss_event_mat = zeros(reps,2);
			
			% loop over repetitions
			for r = 1:reps
				% for repetition r check to see if there were spikes 
				% within event limits
				iei_ind_1 = find((raster{r}>=(neg_lim1(x)*binsize)-binsize)... 
								 & (raster{r}<pos_lim1(x)*binsize));
				if isempty(iei_ind_1)
					iei_ind_1 = 0;
				else
					iei_ind_1 = 1;
				end
									
				iei_ind_2 = find((raster{r}>=(neg_lim2(x)*binsize)-binsize)... 
								 & (raster{r}<pos_lim2(x)*binsize));
				if isempty(iei_ind_2)
					iei_ind_2 = 0;
				else
					iei_ind_2 = 1;
				end
				
				ss_event_mat(r,1) = iei_ind_1;
				ss_event_mat(r,2) = iei_ind_2;
			end % end loop over repetitions
			
			% get percentage of repetitions that have spikes in both 
			% sub-stimulus events
			% ss_event_mat = ss_event_mat(:,1).*ss_event_mat(:,2);
			% ss_events(x) = sum(ss_event_mat)/length(ss_event_mat)*100;
			ss_events(x) = (ss_event_mat(:,1)' * ss_event_mat(:,2)) / reps * 100;			
		end % end loop over sub-stimulus events in cell e
		
		all_ss_events(e).cell_num = displayevents(e).cell_num;
		all_ss_events(e).ss_iei_mat = ss_events;
		% ss_events=[];
    end %if empty ss events
end %for num cells
             
        
if fig
    ss_dist=[];
    for x = 1:ncells
        ss_dist = [ss_dist all_ss_events(x).ss_iei_mat];
    end
    
    [N,X] = hist(ss_dist,100);
    bar(X,N,'b')
    xlabel('Percentage of Overlap in Events')
    ylabel('Count')
    
end
