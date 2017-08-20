function [KernelCorr,KernelPVal] = calckernelcorr(rf1,rf2,stimInfo,Args);
% private kernelcorr function which collapses both RF's into unit vectors
% and calculates the correlation between the two RF's

num_pixels_per_grid_point = stimInfo.data.iniInfo.grid_x_size/stimInfo.data.iniInfo.m_seq_size(1); % Number of pixels in a single M-seq square.
screen_width = stimInfo.data.iniInfo.ScreenWidth;
screen_height = stimInfo.data.iniInfo.ScreenHeight;
if screen_width == 800 & screen_height == 600
    degrees_per_pixel = .05; % Based on the Size of the Monitor, Screen Resolution and Distance of the eyes from the screen.
    pixels_per_degree = 1/.05;
    num_screen_pixels = Args.DegreeDiameter*pixels_per_degree;
    num_kernel_pixels = floor(floor(num_screen_pixels/num_pixels_per_grid_point)/2);
    for i = 1:size(rf1,3);
        rf1_mins(i) = min2(rf1(:,:,i));
        rf1_maxs(i) = max2(rf1(:,:,i));
        rf2_mins(i) = min2(rf2(:,:,i));
        rf2_maxs(i) = max2(rf2(:,:,i));
    end
    % RF1 we find the maximum deviation from mean response of .5
    ranges = [abs(rf1_mins-.5);rf1_maxs-.5];
    [mdr,mdc] = find(ranges == max2(ranges));
    rf_ranges = [rf1_mins;rf1_maxs];
    [rf_center_y,rf_center_x] = find(rf1(:,:,mdc) == rf_ranges(mdr,mdc));
    rf1_center_x = round(mean(rf_center_x));
    rf1_center_y = round(mean(rf_center_y));
    % RF2 we find the maximum deviation from mean response of .5
    ranges = [abs(rf2_mins-.5);rf2_maxs-.5];
    [mdr,mdc] = find(ranges == max2(ranges));
    rf_ranges = [rf2_mins;rf2_maxs];
    [rf_center_y,rf_center_x] = find(rf2(:,:,mdc) == rf_ranges(mdr,mdc));
    rf2_center_x = round(mean(rf_center_x));
    rf2_center_y = round(mean(rf_center_y));
    % Average Center point for both RF's
    rf_center_x = round(mean([rf1_center_x rf2_center_x]));
    rf_center_y = round(mean([rf1_center_y rf2_center_y]));
    %%% Crop out the RF centers for the correlation calculation
    for i = 1:size(rf1,3);
        new_rf1(:,:,i) = rf1(rf_center_y-num_kernel_pixels:rf_center_y+num_kernel_pixels,rf_center_x-num_kernel_pixels:rf_center_x+num_kernel_pixels,i);
        new_rf2(:,:,i) = rf2(rf_center_y-num_kernel_pixels:rf_center_y+num_kernel_pixels,rf_center_x-num_kernel_pixels:rf_center_x+num_kernel_pixels,i);
    end
    
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
    
    % Reshape the kernels into a single unit vector for the correlation
    % calculation
    rf1 = reshape(new_rf1,prod(size(new_rf1)),1);
    rf2 = reshape(new_rf2,prod(size(new_rf2)),1);
    [r,p] = corrcoef(rf1,rf2);
    KernelCorr = r(1,2);
    KernelPVal = p(1,2);
    
    if Args.Plot
        title(['R = ',num2str(KernelCorr) '   P = ',num2str(KernelPVal)])
    end
    
else
    
    KernelCorr=[];
    KernelPVal=[];
    
end