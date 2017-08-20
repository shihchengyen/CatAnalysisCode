function r = plus(p,q,varargin)
%@jointevents/plus Overloaded plus function for Jointevents objects.
%   je = plus(P,Q) combines joint events objects P and Q and returns the
%   JointEvents object R.

% get name of class
classname = mfilename('class');
% check if first input is the right kind of object
if(~isa(p,classname))
    % check if second input is the right kind of object
    if(~isa(q,classname))
        % both inputs are not the right kind of object so create empty
        % object and return it
        r = feval(classname);
    else
        % second input is the right kind of object so return that
        r = q;
    end
else
    if(~isa(q,classname))
        % p is the right kind of object but q is not so just return p
        r = p;
    elseif(isempty(p))
        % p is right object but is empty so return q, which should be
        % right object
        r = q;
    elseif(isempty(q))
        % p are q are both right objects but q is empty while p is not
        % so return p
        r = p;
    else
        % both p and q are the right kind of objects so add them 
        % together
        % assign p to r so that we can be sure we are returning the right
        % object
        r = p;
        r.data.numSets = p.data.numSets + q.data.numSets;
        r.data.setNames = {p.data.setNames{:} q.data.setNames{:}};
        r.data.jointevents.Joint_Overlap_Probability = concatenate(p.data.jointevents.Joint_Overlap_Probability, q.data.jointevents.Joint_Overlap_Probability);
        r.data.jointevents.Expected_Overlap = concatenate(p.data.jointevents.Expected_Overlap, q.data.jointevents.Expected_Overlap);      
        r.data.jointevents.Similarity_Index = concatenate(p.data.jointevents.Similarity_Index, q.data.jointevents.Similarity_Index);
        r.data.jointevents.Joint_Overlap_Durations = concatenate(p.data.jointevents.Joint_Overlap_Durations, q.data.jointevents.Joint_Overlap_Durations);
        r.data.jointevents.Spike_Count_Correlations = concatenate(p.data.jointevents.Spike_Count_Correlations, q.data.jointevents.Spike_Count_Correlations);
        r.data.jointevents.Spike_Count_Significances = concatenate(p.data.jointevents.Spike_Count_Significances, q.data.jointevents.Spike_Count_Significances);  
        r.data.jointevents.Shift_Spike_Count_Correlations = concatenate(p.data.jointevents.Shift_Spike_Count_Correlations, q.data.jointevents.Shift_Spike_Count_Correlations); 
        r.data.jointevents.Shift_Spike_Count_Significances = concatenate(p.data.jointevents.Shift_Spike_Count_Significances, q.data.jointevents.Shift_Spike_Count_Significances);     
        r.data.jointevents.RandomCC_Values = concatenate(p.data.jointevents.RandomCC_Values, q.data.jointevents.RandomCC_Values,'ColumnWise');
        r.data.jointevents.RandomCC_Values_Significance_Value = concatenate(p.data.jointevents.RandomCC_Values_Significance_Value, q.data.jointevents.RandomCC_Values_Significance_Value);   
        r.data.jointevents.Correlation_Coefficent_Value = concatenate(p.data.jointevents.Correlation_Coefficent_Value, q.data.jointevents.Correlation_Coefficent_Value);   
        r.data.jointevents.Correlation_Significance_Value = concatenate(p.data.jointevents.Correlation_Significance_Value, q.data.jointevents.Correlation_Significance_Value);   
        r.data.jointevents.Rasters = concatenate(p.data.jointevents.Rasters , q.data.jointevents.Rasters); 
        r.data.jointevents.Joint_Overlap_Start = concatenate(p.data.jointevents.Joint_Overlap_Start , q.data.jointevents.Joint_Overlap_Start);
        r.data.jointevents.Joint_Overlap_End = concatenate(p.data.jointevents.Joint_Overlap_End , q.data.jointevents.Joint_Overlap_End);        
        r.data.jointevents.NumReps = [p.data.jointevents.NumReps; (p.data.jointevents.NumReps(end) + q.data.jointevents.NumReps)];
        r.data.jointevents.PSTH = concatenate(p.data.jointevents.PSTH , q.data.jointevents.PSTH);
        r.data.jointevents.UpperThresholds = concatenate(p.data.jointevents.UpperThresholds , q.data.jointevents.UpperThresholds);
        r.data.jointevents.LowerThresholds = concatenate(p.data.jointevents.LowerThresholds , q.data.jointevents.LowerThresholds);
        r.data.jointevents.StartEvent = concatenate(p.data.jointevents.StartEvent , q.data.jointevents.StartEvent,'ColumnWise');
        r.data.jointevents.EndEvent = concatenate(p.data.jointevents.EndEvent , q.data.jointevents.EndEvent,'ColumnWise');
        r.nptdata = plus(p.nptdata,q.nptdata);    
    end
end