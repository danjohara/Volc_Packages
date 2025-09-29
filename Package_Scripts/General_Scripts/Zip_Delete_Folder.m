function Zip_Delete_Folder(fol,delFol)
    endChar = fol(end);
    ii = strfind(fol,endChar);
    locFol = fol(ii(end-1)+1:ii(end)-1);
    absZip = [fol(1:end-1),'.zip'];

    if exist(absZip,'file')
        delete(absZip)
    end

    cd(fol); cd ..
    zip([locFol,'.zip'],locFol)
    if delFol
        rmdir(fol,'s')
    end
end