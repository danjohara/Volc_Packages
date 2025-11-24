function Export_Inputs(volcInput,volcPackPrefix)
%%
% Name: Export_Inputs
% Author: Daniel O'Hara
% Data: 05/22/2024 (mm/dd/yyyy)
% Description: Simple function to export algorithm inputs as a text file.
%
% Input:
%   volcInput: Volc package output structure.
%   volcPackPrefix: File prefix to signify which package is calling this
%       function (DrainageVolc, MorVolc, etc.)

if isempty(volcInput.saveResFolder)
    return;
end

fid = fopen([volcInput.saveResFolder,volcInput.figPrefix,volcPackPrefix,'_Input.txt'],'w');
recursiveSave(fid,volcInput,'')
fclose(fid);

end

function recursiveSave(fid,inStruct,dotPref)
    structFields = fieldnames(inStruct);
    
    for i = 1:length(structFields)
        evalc(sprintf('tmpS = inStruct.%s;',structFields{i}));
        if isstruct(tmpS)
            tmpDP = [dotPref,structFields{i},'.'];
            recursiveSave(fid,tmpS,tmpDP);
        else
            if ischar(tmpS)
                fprintf(fid,'%s%s: %s\n',dotPref,structFields{i},tmpS);
            else
                fprintf(fid,'%s%s: %f\n',dotPref,structFields{i},tmpS);
            end
        end
    end
end
            