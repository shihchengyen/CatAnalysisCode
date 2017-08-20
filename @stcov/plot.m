function [obj, varargout] = plot(obj,varargin)
%@stcov/plot Plot function for the stcov object.
%   OBJ = plot(OBJ) creates a stcov plot 

Args = struct('showTitle',1,'xlabel',0,...
    'linkedZoom',0,'showRF',0,'Interp',0);

Args = getOptArgs(varargin,Args,'flags',{'xlabel','linkedZoom','showRF'});

if(isempty(Args.NumericArguments))
    n = 1;
else
    n = Args.NumericArguments{1};
end

%%% Load the eig frames and values
eig_frames = obj.data.eframes(n,:);
eig_values = obj.data.evals(n,:);
st = obj.data.sta(n,:);
%st_f = obj.data.stim_frames(n,:);
for f = 1:size(st,2)
    %sta(:,:,f) = mean(st_f{f},3);
    sta(:,:,f) = st{f};
end
% Calculate the range for the sta
mean_value = mean(sta(:));
if mean_value<1
    sc = max(abs(sta(:)-.5)); %% For the Mseq
    mean_value = .5;
else
    sc = max(abs(sta(:)-mean_value)); %% For the Movies since they are not gamma corrected
    mean_value = mean_value;
end
num_frames = size(eig_frames,2);
count = 0;
if Args.showRF
    for e = 1:num_frames
        e_frames = eig_frames{e};
        e_vals = eig_values{e};
        count=count+1;
        subplot(num_frames,8,count) % 6 eigenframes
        if Args.Interp
            imagesc(interp2(sta(:,:,e)),[(mean_value-sc) (mean_value+sc)])
        else
            imagesc(sta(:,:,e),[(mean_value-sc) (mean_value+sc)])
        end
        axis image
        ylabel(['Frame ' num2str(e)]) 
        if e==1
            title([obj.data.setNames{n}(9:10) 's' obj.data.setNames{n}(16:17) 's' obj.data.setNames{n}(26:27) 'g' obj.data.setNames{n}(37) 'c' obj.data.setNames{n}(47:48)])
        end
        for ee = 1:6
            count=count+1;
            subplot(num_frames,8,count) % 6 eigenframes
            if Args.Interp
                imagesc(interp2(e_frames(:,:,ee)),[min(e_frames(:)) max(e_frames(:))])
            else
                imagesc(e_frames(:,:,ee),[min(e_frames(:)) max(e_frames(:))])
            end
            title(['EIG ',num2str(ee)])
            axis image
            axis off
        end 
        count=count+1;
        subplot(num_frames,8,count) % 6 eigenframes
        plot(1:1:6,e_vals,'-*b');xlim([0 7])
        if e==1
            title('EIG Values')
        end
    end    
else
    for e = 1:num_frames
        e_frames = eig_frames{e};
        e_vals = eig_values{e};
        for ee = 1:6
            count=count+1;
            subplot(num_frames,7,count) % 6 eigenframes
            if Args.Interp
                imagesc(interp2(e_frames(:,:,ee)),[min(e_frames(:)) max(e_frames(:))])
            else
                imagesc(e_frames(:,:,ee),[min(e_frames(:)) max(e_frames(:))])
            end
            title(['EIG ',num2str(ee)])
            axis image
            axis off
        end    
        count=count+1;
        subplot(num_frames,7,count) % 6 eigenframes
        plot(1:1:6,e_vals,'-*b');xlim([0 7])
        if e==1
            title('EIG Values')
        end
    end    
end


