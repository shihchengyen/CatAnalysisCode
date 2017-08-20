s = textread('su3.txt','%s');

for i = 1:length(s)
    cd(s{i})
    [nds,data] = ProcessSession(nptdata,'Cells',ndcells.SessionDirs, ...
        'nptCellCmd','data = adjs2gdf;','DataPlusCmd','data = unique([data; pdata],''Rows'');');

    fname = [strrep(strrep(s{i},'/session','n'),'/site','s') '.gdf'];
    fid = fopen(fname,'wt');
    fprintf(fid,'%i %i\n',round(sortrows(data,2))');
    fclose(fid);
    cd ../../..
end
