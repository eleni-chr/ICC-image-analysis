function ICC_p62puncta
%Loop through subfolders
folders=dir;
idx=[];
for i=1:length(folders)
    if isfolder(folders(i).name)
        idx=[idx;i]; %Ignore items that are not folders
    end
end
folders=folders(idx);

empty=array2table(zeros(0,6));
empty.Properties.VariableNames={'Cell_num','Puncta_num','Area','IntDen','RawIntDen','p62punctaFI'};

punPerCell=array2table(zeros(0,1));
punPerCell.Properties.VariableNames={'PunctaPerCell'};

for j=3:length(folders)
    folname=folders(j).name;
    cd(folname);

    BackgroundIntDenp62=readtable('BackgroundIntDenp62.csv');
    BackgroundIntDenp62.Properties.VariableNames={'ROI_num','Bkg_Area','Bkg_IntDen','Bkg_RawIntDen'};

    p62_Bkg=mean(BackgroundIntDenp62.Bkg_RawIntDen./BackgroundIntDenp62.Bkg_Area);
    
    %Loop through files
    files=dir('p62puncta_*.csv');
    numCells=readtable('p62_CytoplasmIntDen.csv');
    NumOfCells=height(numCells);
    numOfPuncta=zeros(NumOfCells,1);
    cellIDs=(1:NumOfCells)';

    delete 'merged_p62puncta.xlsx' 'p62punctaPerCell.xlsx';

    %Merge info for all cells into a single file
    if ~isempty(files) %ignore folders without the appropriate CSV files
        for m=1:NumOfCells
            writetable(empty,'merged_p62puncta.xlsx','Sheet',num2str(m),'WriteMode','overwritesheet'); %Create an empty Excel sheet for each cell
            writetable(punPerCell,'p62punctaPerCell.xlsx','WriteMode','overwritesheet');
        end
        
        for i=1:length(files)
            data=readtable(files(i).name);
            data.Properties.VariableNames={'Puncta_num','Area','IntDen','RawIntDen'};

            p62punctaFI=(data.RawIntDen-(p62_Bkg.*data.Area))./data.Area;
            for k=1:length(p62punctaFI)
                if p62punctaFI(k)<0
                    p62punctaFI(k)=0; %Adjust negative values
                end
            end
            T2=table(p62punctaFI);

            cellID=extractBetween(files(i).name,'_','.');
            cellNum=str2double(cell2mat(cellID));
            numOfPuncta(cellNum)=data.Puncta_num(end,1); %calculate number of puncta per cell
            Cell_num(1:length(data.Puncta_num))=cellID;
            Cell_num=table(Cell_num');
            data=[Cell_num,data,T2];
            data.Properties.VariableNames{'Var1'}='Cell_num';
            writetable(data,'merged_p62puncta.xlsx','Sheet',cell2mat(cellID),'WriteMode','overwritesheet');
            clear Cell_num
        end
        T=table(cellIDs,numOfPuncta);
        T.Properties.VariableNames={'Cell_num','Total_puncta'};
        writetable(T,'p62punctaPerCell.xlsx','WriteMode','overwritesheet');

        %Combine data into single sheet
        combinedData=[];
        [~,sheets]=xlsfinfo('merged_p62puncta.xlsx');
        for s=1:length(sheets)
            cellData=readtable('merged_p62puncta.xlsx','Sheet',num2str(s));
            if ~isempty(cellData)
                combinedData=[combinedData;cellData];
            end
        end
        writetable(combinedData,'merged_p62puncta.xlsx','Sheet','combined','WriteMode','overwritesheet');
    else %Create empty Excel file for folders without the appropriate CSV files
        for m=1:NumOfCells
            writetable(empty,'merged_p62puncta.xlsx','Sheet',num2str(m),'WriteMode','overwritesheet');
            writetable(punPerCell,'p62punctaPerCell.xlsx','WriteMode','overwritesheet');
        end
    end
    cd ..
end
end