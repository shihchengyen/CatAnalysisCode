function [frames] = ReadPGLFrames(frame_numbers , moviename)
% This function reads the PGL movie files and returns the desired frames
% movienames are Biglebowski, Everest, Animals, Cats, PinkCats, GWN, LSRC
%
% Make sure the movienameroot is correct

movieroot='e:/Natural_Movies/';

for f = 1:length(frame_numbers)

    % Due to the Zero Based Frame numbers of the Presenter code and the
    % One based values in Matlab, we subtract 1. The stiminfo object has
    % the start and end frames corrected for 1 based indexing, hense the
    % subtraction of 1. Scene change files and the frame numbers are 0
    % based.
    frame_number = frame_numbers(f)-1;
    if frame_number < 0
        frame_number = 0;
        fprintf([ 'WARNING: FRAME NUMBER IS LESS THAN ZERO, CHANGING TO ZERO FRAME' '\n'])
    end

    switch(moviename)
        case{'BigLebowski','Biglebowski','biglebowski','pbiglebowski'}
            startstop=[0,3523;3524,7047;7048,10571;10572,13599]; % number of first and last frames in sections of the movie
            movienameroot=[movieroot 'BigLebowski/pbiglebowski00x.pgl']; % Movies must be named as seen above to be read properly.
            for section = 1:size(startstop,1)
                sec = find((frame_number >= startstop(section,1)) & (frame_number <= startstop(section,2)));
                if ~isempty(sec)
                    frame_number = frame_number-startstop(section,1);
                    moviefile=strrep(movienameroot,'x',int2str(section-1));
                    fid=fopen(moviefile, 'r'); %open binary file
                    head=fread(fid, 11, 'char');
                    startframe=fread(fid, 1, 'uint32');
                    lastframe=fread(fid, 1, 'uint32');
                    numframes=fread(fid, 1, 'uint32');
                    fseek(fid,(frame_number)*476*640,0);
                    frame=fread(fid,[640 476*1],'uint8')';
                    fclose(fid);
                end
            end
        case{'Everest','everest'}
            startstop=[0,3494;3495,6989;6990,10484;10485,13979;13980,14999]; % number of first and last frames in sections of the movie
            movienameroot=[movieroot 'Everest/Everest00x.pgl']; % Movies must be named as seen above to be read properly.
            for section = 1:size(startstop,1)
                sec = find((frame_number >= startstop(section,1)) & (frame_number <= startstop(section,2)));
                if ~isempty(sec)
                    frame_number = frame_number-startstop(section,1);
                    moviefile=strrep(movienameroot,'x',int2str(section-1));
                    fid=fopen(moviefile, 'r'); %open binary file
                    head=fread(fid, 11, 'char');
                    startframe=fread(fid, 1, 'uint32');
                    lastframe=fread(fid, 1, 'uint32');
                    numframes=fread(fid, 1, 'uint32');
                    fseek(fid,(frame_number)*480*640,0);
                    frame=fread(fid,[640 480*1],'uint8')';
                    fclose(fid);
                end
            end
        case{'Animals','animals'}
            startstop=[0,3763;3764,7527;7528,11291;11292,14999]; % number of first and last frames in sections of the movie
            movienameroot=[movieroot 'Animals/Animals00x.pgl'];
            for section = 1:size(startstop,1)
                sec = find((frame_number >= startstop(section,1)) & (frame_number <= startstop(section,2)));
                if ~isempty(sec)
                    frame_number = frame_number-startstop(section,1);
                    moviefile=strrep(movienameroot,'x',int2str(section-1));
                    fid=fopen(moviefile, 'r'); %open binary file
                    head=fread(fid, 11, 'char');
                    startframe=fread(fid, 1, 'uint32');
                    lastframe=fread(fid, 1, 'uint32');
                    numframes=fread(fid, 1, 'uint32');
                    fseek(fid,(frame_number)*460*620,0);
                    frame=fread(fid,[620 460*1],'uint8')';
                    fclose(fid);
                end
            end
        case{'Cats','cats'}
            startstop=[0,4556;4557,8999]; % number of first and last frames in sections of the movie
            movienameroot=[movieroot 'Cats/Cats00x.pgl']; % Movies must be named as seen above to be read properly.
            for section = 1:size(startstop,1)
                sec = find((frame_number >= startstop(section,1)) & (frame_number <= startstop(section,2)));
                if ~isempty(sec)
                    frame_number = frame_number-startstop(section,1);
                    moviefile=strrep(movienameroot,'x',int2str(section-1));
                    fid=fopen(moviefile, 'r'); %open binary file
                    head=fread(fid, 11, 'char');
                    startframe=fread(fid, 1, 'uint32');
                    lastframe=fread(fid, 1, 'uint32');
                    numframes=fread(fid, 1, 'uint32');
                    fseek(fid,(frame_number)*380*620,0);
                    frame=fread(fid,[620 380*1],'uint8')';
                    fclose(fid);
                end
            end
        case{'PinkCats','pinkcats'}
            startstop=[0,4556;4557,8999]; % number of first and last frames in sections of the movie
            movienameroot=[movieroot 'PinkCats/PinkCats00x.pgl']; % Movies must be named as seen above to be read properly.
            for section = 1:size(startstop,1)
                sec = find((frame_number >= startstop(section,1)) & (frame_number <= startstop(section,2)));
                if ~isempty(sec)
                    frame_number = frame_number-startstop(section,1);
                    moviefile=strrep(movienameroot,'x',int2str(section-1));
                    fid=fopen(moviefile, 'r'); %open binary file
                    head=fread(fid, 11, 'char');
                    startframe=fread(fid, 1, 'uint32');
                    lastframe=fread(fid, 1, 'uint32');
                    numframes=fread(fid, 1, 'uint32');
                    fseek(fid,(frame_number)*380*620,0);
                    frame=fread(fid,[620 380*1],'uint8')';
                    fclose(fid);
                end
            end
        case{'gwn10'}
            startstop=[0,9999];
            moviefile=[movieroot 'gwn10/gwn10.pgl']; % Movies must be named as seen above to be read properly.
            fid=fopen(moviefile, 'r'); %open binary file
            head=fread(fid, 11, 'char');
            startframe=fread(fid, 1, 'uint32');
            lastframe=fread(fid, 1, 'uint32');
            numframes=fread(fid, 1, 'uint32');
            fseek(fid,(frame_number)*200*200,0);
            frame=fread(fid,[200 200*1],'uint8')';
            fclose(fid);
        case{'gwn15'}
            startstop=[0,4999];
            moviefile=[movieroot 'gwn15/gwn15.pgl']; % Movies must be named as seen above to be read properly.
            fid=fopen(moviefile, 'r'); %open binary file
            head=fread(fid, 11, 'char');
            startframe=fread(fid, 1, 'uint32');
            lastframe=fread(fid, 1, 'uint32');
            numframes=fread(fid, 1, 'uint32');
            fseek(fid,(frame_number)*300*300,0);
            frame=fread(fid,[300 300*1],'uint8')';
            fclose(fid);
        case{'LSRC'}
            startstop=[0,6709;6710,9999]; % number of first and last frames in sections of the movie
            movienameroot=[movieroot 'LSRC/LSRC00x.pgl']; % Movies must be named as seen above to be read properly.
            for section = 1:size(startstop,1)
                sec = find((frame_number >= startstop(section,1)) & (frame_number <= startstop(section,2)));
                if ~isempty(sec)
                    frame_number = frame_number-startstop(section,1);
                    moviefile=strrep(movienameroot,'x',int2str(section-1));
                    fid=fopen(moviefile, 'r'); %open binary file
                    head=fread(fid, 11, 'char');
                    startframe=fread(fid, 1, 'uint32');
                    lastframe=fread(fid, 1, 'uint32');
                    numframes=fread(fid, 1, 'uint32');
                    fseek(fid,(frame_number)*400*400,0);
                    frame=fread(fid,[400 400*1],'uint8')';
                    fclose(fid);
                end
            end

    end %Switch

    frames(:,:,f) = frame;

end