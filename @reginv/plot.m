function [obj, varargout] = plot(obj,varargin)
%@reginv/plot Plot function for the reginv object.
%   OBJ = plot(OBJ) creates a reginv plot 

Args = struct('showRF',0);

Args = getOptArgs(varargin,Args,'flags',{'showRF'});

if(isempty(Args.NumericArguments))
    n = 1;
else
    n = Args.NumericArguments{1};
end

%%% Since the large 4-D matrix contains all the frames, we have to pick out
%%% the frames we want to view.
num_pixels = size(obj.data.reginv(:,:,1,1),2);
num_lamb = size(obj.data.lambda,2);
num_frames = size(obj.data.reginv(:,:,:,1),3);
first = ((num_pixels*n)-num_pixels)+1;
last = num_pixels*n;
for f = 1:num_lamb
    reginv_rf(:,:,:,f) = obj.data.reginv(first:last,:,:,f);
end
count = 1;
for e = 1:num_lamb
    for ee = 1:num_frames
        subplot(num_lamb,num_frames,count)
        % Calculate the range for the sta's
        reginv_mv = mean2(reginv_rf(:,:,:,e)); %% For the Movies since they are not gamma corrected
        sc = abs(reginv_rf(:,:,:,e)-reginv_mv); reginv_sc = max(sc(:)); clear sc
        imagesc(reginv_rf(:,:,ee,e),[(reginv_mv-reginv_sc) (reginv_mv+reginv_sc)]);
        axis image
        if ee == 1
            ylabel([num2str(obj.data.lambda{e})])
        end
        if e == 1
            title(['Frame ',num2str(ee-1)])
            if ee == 1
                title(obj.data.setNames(n))
            end
        end
        count = count + 1;
    end
end



