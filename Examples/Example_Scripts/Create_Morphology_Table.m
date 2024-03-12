%%
% Name: Create_Morphology_Table
% Date: 05/15/2023 (mm/dd/yyyy)
% Author: Daniel O'Hara
% Description: Script to collect DrainageVolc and MorVolc results from 
%   multiple files and save the results to a Matlab and Excel table.
%

%% Setup
packagePath = '.\..\..';
drainageVolc_ResultFolder = '.\..\Example_Results\DrainageVolc\';
morVolc_ResultFolder = '.\..\Example_Results\MorVolc\';

drainageVolc_NameTemplate = 'DrainageVolc_Results.mat';
morVolc_NameTemplate = 'MorVolc_Results.mat';

reliefPercentile = 0.7;

saveFol = '.\..\Example_Results\Compilation\';
saveFilName = 'Example_MatlabTable.mat';
saveExcelName = 'Example_MatlabTable.xlsx';

%% Collect File Lists and Create IDs
addpath(genpath(packagePath))

drainageFiles = ls([drainageVolc_ResultFolder,'*',drainageVolc_NameTemplate]);
morFiles = ls([morVolc_ResultFolder,'*',morVolc_NameTemplate]);

ID_di_mi = {};
for i = 1:size(drainageFiles,1)
    tmp = strtrim(drainageFiles(i,:));
    ID = tmp(1:strfind(tmp,drainageVolc_NameTemplate)-1);
    if strcmp(ID(end),'_')
        ID = ID(1:end-1);
    end

    mi = NaN;
    if contains(morFiles(i,:),ID)
        mi = i;
    else
        for j = 1:size(morFiles,1)
            if contains(morFiles(j,:),ID)
                mi = j;
                break;
            end
        end
    end

    ID_di_mi = [ID_di_mi;{ID,i,mi}];
end

%% Create Table
tmpI = 1;
tableMade = 0;
varNames = {};
varTypes = {};
while ~tableMade
    if isnan(ID_di_mi{tmpI,3})
        tmpI = tmpI + 1;
        continue;
    end

    useDFile = strtrim(drainageFiles(tmpI,:));
    useMFile = strtrim(morFiles(ID_di_mi{tmpI,3},:));

    [dMets,mMets] = Extract_DriainageVolc_MorVolc_Metrics([drainageVolc_ResultFolder,useDFile],[morVolc_ResultFolder,useMFile],reliefPercentile);

    dFields = fieldnames(dMets);
    for i = length(dFields):-1:1
        evalc(sprintf('ss = size(dMets.%s,1);',dFields{i}));
        if ss > 1
            dFields(i) = [];
        end
    end

    mFields = fieldnames(mMets);
    for i = length(mFields):-1:1
        evalc(sprintf('ss = size(mMets.%s,1);',mFields{i}));
        if ss > 1
            mFields(i) = [];
        end
    end

    varNames = ['ID';dFields;mFields];
    varTypes = cell(size(varNames));
    varTypes{1} = 'string';
    for i = 2:length(varNames)
        varTypes{i} = 'double';
    end

    tableMade = 1;
end

Compilation_Table = table('size',[size(drainageFiles,1),length(varNames)],'VariableTypes',varTypes,'VariableNames',varNames);

%% Fill Table
useI = 1;
for i = 1:size(ID_di_mi)
    if isnan(ID_di_mi{i,3})
        warning('Could not find MorVolc file for %s',ID_di_mi{i,1});
        continue;
    end

    useDFile = strtrim(drainageFiles(i,:));
    useMFile = strtrim(morFiles(ID_di_mi{i,3},:));

    [dMets,mMets] = Extract_DriainageVolc_MorVolc_Metrics([drainageVolc_ResultFolder,useDFile],[morVolc_ResultFolder,useMFile],reliefPercentile);
    
    Compilation_Table.ID(useI) = ID_di_mi{useI,1};
    for j = 1:length(dFields)
        try
            evalc(sprintf('Compilation_Table.%s(%d) = dMets.%s;',dFields{j},useI,dFields{j}));
        catch
        end
    end

    for j = 1:length(mFields)
        try
            evalc(sprintf('Compilation_Table.%s(%d) = mMets.%s;',mFields{j},useI,mFields{j}));
        catch
        end
    end

    useI = useI + 1;
end

%% Save Results
save([saveFol,saveFilName],'Compilation_Table');
writetable(Compilation_Table,[saveFol,saveExcelName]);