function showwindowedmovie
%% This function should be called in the session directory to load the
%% stimInfo for the windowed movie sessions or movie sessions
% The Red Points indicate the endpoints of the hand mapped RFs
% The Green Points indicate the size of the window used.

%% Load the stimulus information
st = stiminfo('auto');
if ~strcmp(st.data.iniInfo.type,'Movie')
    fprintf('Not a Windowed Movie Session or Movie Session')
    return
end

win=1;
if ~isfield(st.data.iniInfo,'x_center_small_window_1')
    fprintf('Not a Windowed Movie Session')
    win = 0;
end

if win == 1
    if st.data.iniInfo.enabled_window_1 == 0
        fprintf('Not a Windowed Movie Session')
        if st.data.iniInfo.enabled_window_2 == 1
            win = 0;
        else
            win = 0;
        end
    end
end

if ~isfield(st.data.iniInfo,'ReceptiveField')
    fprintf('No Receptive Field Info')
    return
end

%%%% Get the Screen Dimentions %%%
screen = zeros(st.data.iniInfo.ScreenHeight,st.data.iniInfo.ScreenWidth);
x_offset = (st.data.iniInfo.ScreenWidth - st.data.iniInfo.frame_cols)/2;
y_offset = (st.data.iniInfo.ScreenHeight - st.data.iniInfo.frame_rows)/2;
frame_height = st.data.iniInfo.frame_rows;
frame_width = st.data.iniInfo.frame_cols;
screen(y_offset:frame_height+y_offset-1,x_offset:frame_width+x_offset-1) = 128;
imagesc(screen);colormap gray;hold on; axis image

%%%% Get the RF information %%%%%%%%%
numRF = st.data.iniInfo.ReceptiveField.numRF;
for r = 1:numRF
    x_points = [st.data.iniInfo.ReceptiveField.Points(1:2:end,r); st.data.iniInfo.ReceptiveField.Points(1,r)];
    y_points = [st.data.iniInfo.ReceptiveField.Points(2:2:end,r); st.data.iniInfo.ReceptiveField.Points(2,r)];
    for rr = 1:length(x_points)-1
        plot([x_points(rr) x_points(rr+1)],[y_points(rr) y_points(rr+1)],'r')
    end
end

if win==1
    %%%%%% Plot the Window Movie Information %%%%%%
    if st.data.iniInfo.enabled_window_1 == 1
        plot(st.data.iniInfo.x_center_small_window_1,st.data.iniInfo.y_center_small_window_1,'g*')
        window_radius = st.data.iniInfo.diameter_small_window_1/2;
        plot(st.data.iniInfo.x_center_small_window_1,st.data.iniInfo.y_center_small_window_1+window_radius,'g*')
        plot(st.data.iniInfo.x_center_small_window_1,st.data.iniInfo.y_center_small_window_1-window_radius,'g*')
        plot(st.data.iniInfo.x_center_small_window_1+window_radius,st.data.iniInfo.y_center_small_window_1,'g*')
        plot(st.data.iniInfo.x_center_small_window_1-window_radius,st.data.iniInfo.y_center_small_window_1,'g*')
    end 
    %%%%%% Plot the Window Movie Information %%%%%%
    if st.data.iniInfo.enabled_window_2 == 1
        plot(st.data.iniInfo.x_center_small_window_2,st.data.iniInfo.y_center_small_window_2,'g*')
        window_radius = st.data.iniInfo.diameter_small_window_2/2;
        plot(st.data.iniInfo.x_center_small_window_2+window_radius,st.data.iniInfo.y_center_small_window_2+window_radius,'g*')
        plot(st.data.iniInfo.x_center_small_window_2-window_radius,st.data.iniInfo.y_center_small_window_2+window_radius,'g*')
        plot(st.data.iniInfo.x_center_small_window_2+window_radius,st.data.iniInfo.y_center_small_window_2-window_radius,'g*')
        plot(st.data.iniInfo.x_center_small_window_2-window_radius,st.data.iniInfo.y_center_small_window_2-window_radius,'g*')
    end
end