% rf_info.m - get rf coords in mseq array
%
% function [rmin rmax cmin cmax] = rf_info(session_name,group,cell)

function [rmin, rmax, cmin, cmax, f_nums, max_frame] = rf_info(session_name);

switch session_name
    case 't107'
        rmin=25;
        rmax=31;
        cmin=25;
        cmax=32;
        f_nums = [2 6];
        max_frame = 2;
    case 't205'
        rmin=25;
        rmax=33;
        cmin=29;
        cmax=38;
        f_nums = [2 6];
        max_frame = 2;
    case 't214'
        rmin=28;
        rmax=35;
        cmin=26;
        cmax=33; 
        f_nums = [2 6];
        max_frame = 2;
    case 'a225'
        rmin=27;
        rmax=36;
        cmin=29;
        cmax=39;
        f_nums = [2 5];
        max_frame = 3;
    case 'a231'
        rmin=7;
        rmax=13;
        cmin=6;
        cmax=13;
        f_nums = [2 6];
        max_frame = 2;
    case 'a412'
        rmin=22;
        rmax=43;
        cmin=25;
        cmax=46;   
        f_nums = [2 5];
        max_frame = 2;
    case 'a416'
        rmin=24;
        rmax=39;
        cmin=21;
        cmax=38;
        f_nums = [2 4];
        max_frame = 2;
    case 'a426'        
        rmin=21;
        rmax=37;
        cmin=24;
        cmax=37;   
        f_nums = [1 5];
        max_frame = 2;
end

