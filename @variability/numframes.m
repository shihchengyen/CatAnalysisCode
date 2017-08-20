    aframes = Args.NumFrames;
    % check division of numframes into sub-bins
    sbinsperframe = framebins/aframes;
    if(rem(sbinsperframe,1)==0)
        % we can just do one histogram and reuse the values
	    dafv = diff(afv);
        binsizes = dafv/sbinsperframe;
		mat1 = tril(ones(sbinsperframe));
		mat2 = [afv(1:(afvl-1))'; repmat(binsizes',(sbinsperframe-1),1)];
		blimits = mat1 * mat2;
		binlimits = [blimits(:); afv(afvl)];
		sc = histcie(data.cell_info.adjusted_spiketrain, ...
			binlimits,'DropLast');
        if(aframes>1)
			% reshape into framebins x (aframes * nreps)
            nfullframes = nframes - aframes + 1;
            ntotalframes = nfullframes * repetitions;
            i1 = reshape(0:(nframes*repetitions-1),nframes,[]);
            i2 = i1(1:nfullframes,:) * sbinsperframe + 1;
            i3 = reshape(i2,1,[]);
            i4 = tril(ones(framebins));
            i5 = [i3; ones(framebins-1,ntotalframes)];
            i6 = i4 * i5;
            sc1 = sc(i6);
        else
    		sc1 = reshape(sc,framebins,[]);
        end
    else
        % we will have to do multiple histograms

        
        if(aframes>1)
            
        else
        end
        % get integral multiple of aframes
        afvl2 = floor(afvl/aframes);
        afvl2e = afvl2*aframes;
        afv2 = reshape(afv(1:afvl2e),aframes,afvl2);
        mat2r1 = afv2(1,:);
        dafv = [diff(mat2r1) 

        dafv = diff(afv);
        mat2r1 = afv(1:(afvl-1))';
        
		binsizes = dafv/framebins;
		mat1 = tril(ones(framebins));
		mat2 = [mat2r1; repmat(binsizes',(framebins-1),1)];
    end
