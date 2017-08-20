function [KernelCorr,KernelPVal] = brunokernelcorr(rf1,rf2,session_name,Args)
% This function does the same calculation but the cropping is done by a
% function of known edge points.

[rmin, rmax, cmin, cmax, f_nums, max_frame]=rf_info(session_name);
f_nums = f_nums(1):f_nums(2); % Correlate to the Best Kernels over time
%f_nums = max_frame; % Correlate to the Maximum Kernel Frame only.
for i = 1:length(f_nums)
    if Args.Threshold
        high_thresh = prctile(rf1(:),90);
        low_thresh = prctile(rf1(:),10);    
        n_rf1 = rf1(rmin:rmax,cmin:cmax,f_nums(i));
        index = find(n_rf1(:,:,i)<high_thresh & n_rf1(:,:,i)>low_thresh);
        n_rf1(index) = 0;
        new_rf1(:,:,i) = n_rf1;
        high_thresh = prctile(rf2(:),90);
        low_thresh = prctile(rf2(:),10);    
        n_rf2 = rf2(rmin:rmax,cmin:cmax,f_nums(i));
        index = find(n_rf2(:,:,i)<high_thresh & n_rf2(:,:,i)>low_thresh);
        n_rf2(index) = 0;
        new_rf2(:,:,i) = n_rf2;
    else
        new_rf1(:,:,i) = rf1(rmin:rmax,cmin:cmax,f_nums(i));
        new_rf2(:,:,i) = rf2(rmin:rmax,cmin:cmax,f_nums(i));
    end
end
% Reshape the kernels into a single unit vector for the correlation
% calculation
if Args.Plot
    minv_1 = .5-min(rf1(:));
    maxv_1 = max(rf1(:))-.5;
    maxv_1 = max([minv_1 maxv_1]);
    minv_1 = .5-maxv_1;
    maxv_1 = .5+maxv_1;    
    minv_2 = .5-min(rf2(:));
    maxv_2 = max(rf2(:))-.5;
    maxv_2 = max([minv_2 maxv_2]);
    minv_2 = .5-maxv_2;
    maxv_2 = .5+maxv_2; 
    num_plots = size(new_rf1,3);
    scrsz = get(0,'ScreenSize'); % This is so that the figure will nicely fit into the powerpoint presentation
    figure('Position',[200 200 scrsz(3)*.715 scrsz(4)/1.5])
    for f = 1:num_plots
        subplot(2,num_plots,f)
        imagesc(new_rf1(:,:,f),[minv_1 maxv_1]); axis image;
        if f == 1
            title(Args.ClusterDirs{1})
        end
    end
    for f = 1:num_plots
        subplot(2,num_plots,f+num_plots)
        imagesc(new_rf2(:,:,f),[minv_2 maxv_2]); axis image;
        if f == 1
            title(Args.ClusterDirs{2})
        end
    end
end

rf1 = reshape(new_rf1,prod(size(new_rf1)),1);
rf2 = reshape(new_rf2,prod(size(new_rf2)),1);
[r,p] = corrcoef(rf1,rf2);
KernelCorr = r(1,2);
KernelPVal = p(1,2);

if Args.Plot
    title(['R = ',num2str(KernelCorr) '   P = ',num2str(KernelPVal)])
end