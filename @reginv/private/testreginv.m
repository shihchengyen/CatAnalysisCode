function [R,S] = testreginv(Args,Threshold,movie_frame_numbers,moviename,center_point,screen_height,screen_width,y_off,x_off,frame_height,frame_width,num_frames);
    
% Get the gabor function to be convolved with the movie
gv=gaborphase(20,15,15,Args.NumPixels); 
sum_gv = sum(sum(gv.*gv));

if Args.Reps
    movie_frame_numbers = movie_frame_numbers(1:Args.Reps*num_frames);
end

% Convolve the gabor function with the Movie Frames
num_frames_back = (1:Args.NumFramesBack)-1;
pix = Args.NumPixels/2;
for ii = 1:Args.NumFramesBack
    movie_frame_numbers = movie_frame_numbers-num_frames_back(ii);
    SS = uint8(zeros(length(movie_frame_numbers),Args.NumPixels^2)); %% Put into Memory the S Matrix
    for ff = 1:length(movie_frame_numbers)
        movie_frame = ReadPGLFrames(movie_frame_numbers(ff),moviename);
        screen_image = zeros(screen_height,screen_width);
        screen_image(y_off:frame_height+y_off-1,x_off:frame_width+x_off-1) = movie_frame;
        screen_image = (screen_image(center_point(1)-(pix-1):center_point(1)+pix,center_point(2)-(pix-1):center_point(2)+pix));
        SS(ff,:) = reshape(screen_image,1,Args.NumPixels^2);
        if ii == 1;
            if strcmp(Args.TestType,'Random')
                r = rand;
                if r >= .5
                    R(ff) = 1;
                else
                    R(ff) = 0;
                end
            elseif strcmp(Args.TestType,'Convolve')            
                c = corr2(screen_image,gv);
                if c >= (Threshold/100) 
                    R(ff) = 1;
                else
                    R(ff) = 0;
                end                 
            else
                conv_image = screen_image.*gv;                
                if sum(conv_image(:))>=((sum_gv)*(Threshold/100))
                    R(ff) = 1;
                else
                    R(ff) = 0;
                end                 
            end
        end 
        S{ii} = SS;
    end
    fprintf(['Frames Back ',num2str(ii) '\n'])
end

if strcmp(Args.TestType,'Convolve')
    ind = find(R == 1);
    R = R(ind);
    for ii = 1:Args.NumFramesBack
        ss{ii} = S{ii}(ind,:);        
    end
    S = ss;
end
if strcmp(Args.TestType,'Random')
    %rep = R(1:num_frames); R = [];
    %for r = 1:Args.Reps
    %    R = [R rep];
    %end    
end
R = R';