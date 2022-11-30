function ICC_merge_data

sampleNTG={'163B2','163B3','163B4','167B2','167B3','167B4','172A1','172A2','172A3'};
sampleLOA={'165B1','171A1','171A2','171A3','171C1','171C2','171C3'};
sampleTGW={'180A1','180B1','180A2','180B2','180A5','180B5'};
sampleTGM={'176A1','176A2','176A3','177A1','177A2','177A3','182A3','182B3'};
sampleLWT={'178C3','178D3','179A1','179A2'};
sampleLMT={'164B3','164B4','164B5','182A1','182B1','182A2','182B2'};

mouseNTG1={'163B2','163B3','163B4'};
mouseNTG2={'167B2','167B3','167B4'};
mouseNTG3={'172A1','172A2','172A3'};
mouseLOA1={'165B1'};
mouseLOA2={'171A1','171A2','171A3'};
mouseLOA3={'171C1','171C2','171C3'};
mouseTGW1={'180A1','180B1'};
mouseTGW2={'180A2','180B2'};
mouseTGW3={'180A5','180B5'};
mouseTGM1={'176A1','176A2','176A3'};
mouseTGM2={'177A1','177A2','177A3'};
mouseTGM3={'182A3','182B3'};
mouseLWT1={'178C3','178D3'};
mouseLWT2={'179A1'};
mouseLWT3={'179A2'};
mouseLMT1={'164B3','164B4','164B5'};
mouseLMT2={'182A1','182B1'};
mouseLMT3={'182A2','182B2'};

parentdir=cd;

%% Analyse images in Batch 1
cd("Batch 1")
%Loop through subfolders
folders=dir;
for j=3:length(folders)
    if ~isfolder(folders(j).name)
        continue; %Ignore items that are not folders
    end
    
    folname=folders(j).name;
    if strcmp(folname,'Slide207prox_172A2_10')
        continue;
    end
    cd(folname);
    
    %Load data
    BackgroundIntDenTDP=readtable('BackgroundIntDenTDP.csv');
    BackgroundIntDenTDP.Properties.VariableNames={'ROI_num','Bkg_Area','Bkg_IntDen','Bkg_RawIntDen'};
    
    BackgroundIntDenp62=readtable('BackgroundIntDenp62.csv');
    BackgroundIntDenp62.Properties.VariableNames={'ROI_num','Bkg_Area','Bkg_IntDen','Bkg_RawIntDen'};    

    TDP_NucleusIntDen=readtable('TDP_NucleusIntDen.csv');
    TDP_NucleusIntDen.Properties.VariableNames={'Cell_num','Area_Nuc','TDP_IntDen_Nuc','TDP_RawIntDen_Nuc'};
    
    TDP_CytoplasmIntDen=readtable('TDP_CytoplasmIntDen.csv');
    TDP_CytoplasmIntDen.Properties.VariableNames={'Cell_num','Area_Cyt','TDP_IntDen_Cyt','TDP_RawIntDen_Cyt'};
    
    p62cytoplasmIntDen=readtable('p62_CytoplasmIntDen.csv');
    p62cytoplasmIntDen.Properties.VariableNames={'Cell_num','p62_Area_Cyt','p62_IntDen_Cyt','p62_RawIntDen_Cyt'};

    p62punctaPerCell=readtable('p62punctaPerCell.xlsx');
    if ~isempty(p62punctaPerCell)
        p62punctaPerCell.Properties.VariableNames={'Cell_num','Total_p62_puncta'};
    else
        p62punctaPerCell=array2table(nan(height(TDP_NucleusIntDen),1));
        p62punctaPerCell.Properties.VariableNames={'Total_p62_puncta'};
    end
    
    %Concatenate data into table
    Tbl=[TDP_NucleusIntDen,TDP_CytoplasmIntDen(:,2:end),p62cytoplasmIntDen(:,2:end)];
    Tbl_final=[TDP_NucleusIntDen(:,1:2),TDP_CytoplasmIntDen(:,2)];

    %Calculations for TDP-43
    TDP_Bkg=mean(BackgroundIntDenTDP.Bkg_RawIntDen./BackgroundIntDenTDP.Bkg_Area); %background fluorescence intensity in TDP-43 channel
    TDP_N=(Tbl.TDP_RawIntDen_Nuc-(TDP_Bkg.*Tbl.Area_Nuc))./Tbl.Area_Nuc; %TDP-43 fluorescence intensity in the nucleus
    TDP_C=(Tbl.TDP_RawIntDen_Cyt-(TDP_Bkg.*Tbl.Area_Cyt))./Tbl.Area_Cyt; %TDP-43 fluorescence intensity in the cytoplasm

    %Adjust negative values (i.e., when there is less signal in the ROI than in the background)
    for i=1:length(TDP_N)
        if TDP_N(i)<0
            TDP_N(i)=0;
        end
        if TDP_C(i)<0
            TDP_C(i)=0;
        end
    end

    TDP_T=TDP_N+TDP_C; %TDP-43 fluorescence intensity in soma
    TDP_pN=TDP_N./TDP_T.*100; %percentage of TDP-43 in nucleus
    TDP_pC=TDP_C./TDP_T.*100; %percentage of TDP-43 in cytoplasm
    
    %Calculate if the cell contains more than 50% cytoplasmic TDP-43
    MajorCytTDP=zeros(length(TDP_pC),1);
    for i=1:length(TDP_pC)
        if TDP_pC(i)>50
            MajorCytTDP(i)=1; %majority of TDP-43 is in cytoplasm
        else
            MajorCytTDP(i)=0; %majority of TDP-43 is in nucleus
        end
    end
    
    %Calculations for p62
    p62_Bkg=mean(BackgroundIntDenp62.Bkg_RawIntDen./BackgroundIntDenp62.Bkg_Area); %background fluorescence intensity in p62 channel
    p62_C=(Tbl.p62_RawIntDen_Cyt-(p62_Bkg.*Tbl.p62_Area_Cyt))./Tbl.p62_Area_Cyt; %p62 fluorescence intensity in cytoplasm

    %Adjust negative values (i.e., when there is less signal in the ROI than in the background)
    for i=1:length(p62_C)
        if p62_C(i)<0
            p62_C(i)=0;
        end
    end

    %Normalise puncta number to cytoplasm area
    p62punPerCell=p62punctaPerCell.Total_p62_puncta./Tbl.Area_Cyt;

    %Concatenate calculation results into table
    Tbl2=array2table([TDP_N,TDP_C,TDP_T,TDP_pN,TDP_pC,MajorCytTDP,p62_C,p62punPerCell],...
        'VariableNames',{'TDP_Nuc_FI','TDP_Cyt_FI','TDP_Soma_FI','TDP_percent_Nuc','TDP_percent_Cyt','MajorCytTDP','p62_Cyt_FI','Total_p62_puncta'});

    %Concatenate raw data and calculations into table
    Data=[Tbl_final,Tbl2];
    
    %Extract identifier info and append to table
    seq=extractAfter(folname,'_');
    ID=extractBefore(seq,'_');
    sampleID=char(zeros(length(TDP_N),length(ID)));
    mouseID=char(zeros(length(TDP_N),4));
    genotype=char(zeros(length(TDP_N),24));
    for i=1:length(TDP_N)
        sampleID(i,:)=ID;
        if ismember(sampleID(i,:),sampleNTG)
            genotype(i,:)='Dync1h1(+/+),TARDBP(-/-)';
            if ismember(ID,mouseNTG1)
                mouseID(i,:)='NTG1';
            elseif ismember(ID,mouseNTG2)
                mouseID(i,:)='NTG2';
            elseif ismember(ID,mouseNTG3)
                mouseID(i,:)='NTG3';
            end
        elseif ismember(sampleID(i,:),sampleLOA)
            genotype(i,:)='Dync1h1(+/L),TARDBP(-/-)';
            if ismember(ID,mouseLOA1)
                mouseID(i,:)='LOA1';
            elseif ismember(ID,mouseLOA2)
                mouseID(i,:)='LOA2';
            elseif ismember(ID,mouseLOA3)
                mouseID(i,:)='LOA3';
            end
        end
    end
    Data.sampleID=sampleID;
    Data.genotype=genotype;
    Data.mouseID=mouseID;
    
    image_num=extractAfter(seq,'_');
    if length(image_num)==1
        image_num=strcat('0',image_num);
    end
    im_num=char(zeros(length(TDP_N),2));
    for i=1:length(TDP_N)
        im_num(i,:)=image_num;
    end
    Data.image_num=im_num;
    
    writetable(Data,'mergedData.xlsx','WriteMode','replacefile');
    cd ..
end
cd(parentdir)

%% Analyse images in Batch 2
cd("Batch 2")
%Loop through subfolders
folders=dir;
for j=3:length(folders)
    if ~isfolder(folders(j).name)
        continue; %Ignore items that are not folders
    end
    
    folname=folders(j).name;
    cd(folname);
    
    %Load data
    BackgroundIntDenTDP=readtable('BackgroundIntDenTDP.csv');
    BackgroundIntDenTDP.Properties.VariableNames={'ROI_num','Bkg_Area','Bkg_IntDen','Bkg_RawIntDen'};
    
    BackgroundIntDenp62=readtable('BackgroundIntDenp62.csv');
    BackgroundIntDenp62.Properties.VariableNames={'ROI_num','Bkg_Area','Bkg_IntDen','Bkg_RawIntDen'};    

    BackgroundIntDenGFP=readtable('BackgroundIntDenGFP.csv');
    BackgroundIntDenGFP.Properties.VariableNames={'ROI_num','Bkg_Area','Bkg_IntDen','Bkg_RawIntDen'};

    TDP_NucleusIntDen=readtable('TDP_NucleusIntDen.csv');
    TDP_NucleusIntDen.Properties.VariableNames={'Cell_num','Area_Nuc','TDP_IntDen_Nuc','TDP_RawIntDen_Nuc'};
    
    TDP_CytoplasmIntDen=readtable('TDP_CytoplasmIntDen.csv');
    TDP_CytoplasmIntDen.Properties.VariableNames={'Cell_num','Area_Cyt','TDP_IntDen_Cyt','TDP_RawIntDen_Cyt'};
    
    p62cytoplasmIntDen=readtable('p62_CytoplasmIntDen.csv');
    p62cytoplasmIntDen.Properties.VariableNames={'Cell_num','p62_Area_Cyt','p62_IntDen_Cyt','p62_RawIntDen_Cyt'};

    GFP_NucleusIntDen=readtable('TDP_NucleusIntDen.csv');
    GFP_NucleusIntDen.Properties.VariableNames={'Cell_num','Area_Nuc','GFP_IntDen_Nuc','GFP_RawIntDen_Nuc'};
    
    GFP_CytoplasmIntDen=readtable('TDP_CytoplasmIntDen.csv');
    GFP_CytoplasmIntDen.Properties.VariableNames={'Cell_num','Area_Cyt','GFP_IntDen_Cyt','GFP_RawIntDen_Cyt'};
    
    p62punctaPerCell=readtable('p62punctaPerCell.xlsx');
    if ~isempty(p62punctaPerCell)
        p62punctaPerCell.Properties.VariableNames={'Cell_num','Total_p62_puncta'};
    else
        p62punctaPerCell=array2table(nan(height(TDP_NucleusIntDen),1));
        p62punctaPerCell.Properties.VariableNames={'Total_p62_puncta'};
    end

    %Concatenate data into table
    Tbl=[TDP_NucleusIntDen,TDP_CytoplasmIntDen(:,2:end),p62cytoplasmIntDen(:,2:end),GFP_NucleusIntDen(:,3:end),GFP_CytoplasmIntDen(:,3:end)];
    Tbl_final=[TDP_NucleusIntDen(:,1:2),TDP_CytoplasmIntDen(:,2)];

    %Calculations for TDP-43
    TDP_Bkg=mean(BackgroundIntDenTDP.Bkg_RawIntDen./BackgroundIntDenTDP.Bkg_Area); %background fluorescence intensity in TDP-43 channel
    TDP_N=(Tbl.TDP_RawIntDen_Nuc-(TDP_Bkg.*Tbl.Area_Nuc))./Tbl.Area_Nuc; %TDP-43 fluorescence intensity in the nucleus
    TDP_C=(Tbl.TDP_RawIntDen_Cyt-(TDP_Bkg.*Tbl.Area_Cyt))./Tbl.Area_Cyt; %TDP-43 fluorescence intensity in the cytoplasm

    %Adjust negative values (i.e., when there is less signal in the ROI than in the background)
    for i=1:length(TDP_N)
        if TDP_N(i)<0
            TDP_N(i)=0;
        end
        if TDP_C(i)<0
            TDP_C(i)=0;
        end
    end

    TDP_T=TDP_N+TDP_C; %TDP-43 fluorescence intensity in soma
    TDP_pN=TDP_N./TDP_T.*100; %percentage of TDP-43 in nucleus
    TDP_pC=TDP_C./TDP_T.*100; %percentage of TDP-43 in cytoplasm
    
    %Calculate if the cell contains more than 50% cytoplasmic TDP-43
    MajorCytTDP=zeros(length(TDP_pC),1);
    for i=1:length(TDP_pC)
        if TDP_pC(i)>50
            MajorCytTDP(i)=1; %majority of TDP-43 is in cytoplasm
        else
            MajorCytTDP(i)=0; %majority of TDP-43 is in nucleus
        end
    end

    %Calculations for p62
    p62_Bkg=mean(BackgroundIntDenp62.Bkg_RawIntDen./BackgroundIntDenp62.Bkg_Area); %background fluorescence intensity in p62 channel
    p62_C=(Tbl.p62_RawIntDen_Cyt-(p62_Bkg.*Tbl.p62_Area_Cyt))./Tbl.p62_Area_Cyt; %p62 fluorescence intensity in cytoplasm

    %Adjust negative values (i.e., when there is less signal in the ROI than in the background)
    for i=1:length(p62_C)
        if p62_C(i)<0
            p62_C(i)=0;
        end
    end

    %Normalise puncta number to cytoplasm area
    p62punPerCell=p62punctaPerCell.Total_p62_puncta./Tbl.Area_Cyt;

    %Calculations for GFP
    GFP_Bkg=mean(BackgroundIntDenGFP.Bkg_RawIntDen./BackgroundIntDenGFP.Bkg_Area); %background fluorescence intensity in GFP channel
    GFP_N=(Tbl.GFP_RawIntDen_Nuc-(GFP_Bkg.*Tbl.Area_Nuc))./Tbl.Area_Nuc; %GFP fluorescence intensity in the nucleus
    GFP_C=(Tbl.GFP_RawIntDen_Cyt-(GFP_Bkg.*Tbl.Area_Cyt))./Tbl.Area_Cyt; %GFP fluorescence intensity in the cytoplasm

    %Adjust negative values (i.e., when there is less signal in the ROI than in the background)
    for i=1:length(GFP_N)
        if GFP_N(i)<0
            GFP_N(i)=0;
        end
        if GFP_C(i)<0
            GFP_C(i)=0;
        end
    end

    GFP_T=GFP_N+GFP_C; %GFP fluorescence intensity in soma
    GFP_pN=GFP_N./GFP_T.*100; %percentage of GFP in nucleus
    GFP_pC=GFP_C./GFP_T.*100; %percentage of GFP in cytoplasm
    
    %Calculate if the cell contains more than 50% cytoplasmic GFP
    MajorCytGFP=zeros(length(GFP_pC),1);
    for i=1:length(GFP_pC)
        if GFP_pC(i)>50
            MajorCytGFP(i)=1; %majority of GFP is in cytoplasm
        else
            MajorCytGFP(i)=0; %majority of GFP is in nucleus
        end
    end

    %Concatenate calculation results into table
    Tbl2=array2table([TDP_N,TDP_C,TDP_T,TDP_pN,TDP_pC,MajorCytTDP,p62_C,GFP_N,GFP_C,GFP_T,GFP_pN,GFP_pC,MajorCytGFP,p62punPerCell],...
        'VariableNames',{'TDP_Nuc_FI','TDP_Cyt_FI','TDP_Soma_FI','TDP_percent_Nuc','TDP_percent_Cyt','MajorCytTDP',...
        'p62_Cyt_FI','GFP_Nuc_FI','GFP_Cyt_FI','GFP_Soma_FI','GFP_percent_Nuc','GFP_percent_Cyt','MajorCytGFP','Total_p62_puncta'});

    %Concatenate raw data and calculations into table
    Data=[Tbl_final,Tbl2];
    
    %Extract identifier info and append to table
    seq=extractAfter(folname,'_');
    ID=extractBefore(seq,'_');
    sampleID=char(zeros(length(TDP_N),length(ID)));
    mouseID=char(zeros(length(TDP_N),4));
    genotype=char(zeros(length(TDP_N),24));
    for i=1:length(TDP_N)
        sampleID(i,:)=ID;
        if ismember(sampleID(i,:),sampleTGW)
            genotype(i,:)='Dync1h1(+/+),TARDBP(-/+)';
            if ismember(ID,mouseTGW1)
                mouseID(i,:)='TGW1';
            elseif ismember(ID,mouseTGW2)
                mouseID(i,:)='TGW2';
            elseif ismember(ID,mouseTGW3)
                mouseID(i,:)='TGW3';
            end
        elseif ismember(sampleID(i,:),sampleTGM)
            genotype(i,:)='Dync1h1(+/+),TARDBP(-/M)';
            if ismember(ID,mouseTGM1)
                mouseID(i,:)='TGM1';
            elseif ismember(ID,mouseTGM2)
                mouseID(i,:)='TGM2';
            elseif ismember(ID,mouseTGM3)
                mouseID(i,:)='TGM3';
            end
        elseif ismember(sampleID(i,:),sampleLWT)
            genotype(i,:)='Dync1h1(+/L),TARDBP(-/+)';
            if ismember(ID,mouseLWT1)
                mouseID(i,:)='LWT1';
            elseif ismember(ID,mouseLWT2)
                mouseID(i,:)='LWT2';
            elseif ismember(ID,mouseLWT3)
                mouseID(i,:)='LWT3';
            end
        elseif ismember(sampleID(i,:),sampleLMT)
            genotype(i,:)='Dync1h1(+/L),TARDBP(-/M)';
            if ismember(ID,mouseLMT1)
                mouseID(i,:)='LMT1';
            elseif ismember(ID,mouseLMT2)
                mouseID(i,:)='LMT2';
            elseif ismember(ID,mouseLMT3)
                mouseID(i,:)='LMT3';
            end
        end
    end
    Data.sampleID=sampleID;
    Data.genotype=genotype;
    Data.mouseID=mouseID;
    
    image_num=extractAfter(seq,'_');
    if length(image_num)==1
        image_num=strcat('0',image_num);
    end
    im_num=char(zeros(length(TDP_N),2));
    for i=1:length(TDP_N)
        im_num(i,:)=image_num;
    end
    Data.image_num=im_num;
    
    writetable(Data,'mergedData.xlsx','WriteMode','replacefile');
    cd ..
end
cd(parentdir)

%% Create master file
cd("Batch 1")
%Loop through subfolders
folders=dir;
idx=[];
for i=1:length(folders)
    if isfolder(folders(i).name)
        idx=[idx;i]; %Ignore items that are not folders
    end
end
folders=folders(idx);

T=[];
for j=3:length(folders)
    folname=folders(j).name;
    if strcmp(folname,'Slide207prox_172A2_10')
        continue;
    end
    cd(folname);
    
    if isfile('mergedData.xlsx')
        data=readtable('mergedData.xlsx');
        T=[T;data];
    end
    cd ..
end
cd(parentdir)
writetable(T,'masterfile.xlsx','Sheet','Batch 1','WriteMode','overwritesheet');

% Add additional sheets to the Excel file with separate genotypes
NTG=T(strcmp(T.genotype,'Dync1h1(+/+),TARDBP(-/-)'),:);
LOA=T(strcmp(T.genotype,'Dync1h1(+/L),TARDBP(-/-)'),:);

writetable(NTG,'masterfile.xlsx','Sheet','NTG','WriteMode','overwritesheet');
writetable(LOA,'masterfile.xlsx','Sheet','LOA','WriteMode','overwritesheet');

% Additional calculations per animal
T_NTG=[];
T_LOA=[];

for i=1:3
    str=strcat('NTG',num2str(i));
    genotype={'Dync1h1(+/+),TARDBP(-/-)'};
    lookup=NTG(strcmp(NTG.mouseID,str),:);
    percent_cells_major_cyt_TDP=num2cell(sum(lookup.MajorCytTDP)/height(lookup)*100); %percentage of cells with majority cytoplasmic TDP-43
    totalp62=height(lookup)-sum(isnan(lookup.Total_p62_puncta)); %total number of cells (without the excluded cells)
    percent_cells_p62_puncta=num2cell((length(find(lookup.Total_p62_puncta))-sum(isnan(lookup.Total_p62_puncta)))/totalp62*100); %percentage of cells with p62 puncta
    T_NTG=[T_NTG;table({str},genotype,percent_cells_major_cyt_TDP,percent_cells_p62_puncta)];
end
T_NTG.Properties.VariableNames={'mouseID','genotype','percent_cells_MajorCytTDP','percent_cells_p62_puncta'};

for i=1:3
    str=strcat('LOA',num2str(i));
    genotype={'Dync1h1(+/L),TARDBP(-/-)'};
    lookup=LOA(strcmp(LOA.mouseID,str),:);
    percent_cells_major_cyt_TDP=num2cell(sum(lookup.MajorCytTDP)/height(lookup)*100); %percentage of cells with majority cytoplasmic TDP-43
    totalp62=height(lookup)-sum(isnan(lookup.Total_p62_puncta)); %total number of cells (without the excluded cells)
    percent_cells_p62_puncta=num2cell((length(find(lookup.Total_p62_puncta))-sum(isnan(lookup.Total_p62_puncta)))/totalp62*100); %percentage of cells with p62 puncta
    T_LOA=[T_LOA;table({str},genotype,percent_cells_major_cyt_TDP,percent_cells_p62_puncta)];
end
T_LOA.Properties.VariableNames={'mouseID','genotype','percent_cells_MajorCytTDP','percent_cells_p62_puncta'};

T2=[T_NTG;T_LOA];
writetable(T2,'masterfile.xlsx','Sheet','Additional data Batch 1','WriteMode','overwritesheet');

cd("Batch 2")
%Loop through subfolders
folders=dir;
idx=[];
for i=1:length(folders)
    if isfolder(folders(i).name)
        idx=[idx;i]; %Ignore items that are not folders
    end
end
folders=folders(idx);

T=[];
for j=3:length(folders)
    folname=folders(j).name;
    cd(folname);
    
    if isfile('mergedData.xlsx')
        data=readtable('mergedData.xlsx');
        T=[T;data];
    end
    cd ..
end
cd(parentdir)
writetable(T,'masterfile.xlsx','Sheet','Batch 2','WriteMode','overwritesheet');

% Add additional sheets to the Excel file with separate genotypes
TGW=T(strcmp(T.genotype,'Dync1h1(+/+),TARDBP(-/+)'),:);
TGM=T(strcmp(T.genotype,'Dync1h1(+/+),TARDBP(-/M)'),:);
LWT=T(strcmp(T.genotype,'Dync1h1(+/L),TARDBP(-/+)'),:);
LMT=T(strcmp(T.genotype,'Dync1h1(+/L),TARDBP(-/M)'),:);

writetable(TGW,'masterfile.xlsx','Sheet','WT-Tg','WriteMode','overwritesheet');
writetable(TGM,'masterfile.xlsx','Sheet','MUT-Tg','WriteMode','overwritesheet');
writetable(LWT,'masterfile.xlsx','Sheet','LOA-WT-Tg','WriteMode','overwritesheet');
writetable(LMT,'masterfile.xlsx','Sheet','LOA-MUT-Tg','WriteMode','overwritesheet');

% Additional calculations per animal
T_TGW=[];
T_TGM=[];
T_LWT=[];
T_LMT=[];

for i=1:3
    str=strcat('TGW',num2str(i));
    genotype={'Dync1h1(+/+),TARDBP(-/+)'};
    lookup=TGW(strcmp(TGW.mouseID,str),:);
    percent_cells_major_cyt_TDP=num2cell(sum(lookup.MajorCytTDP)/height(lookup)*100); %percentage of cells with majority cytoplasmic TDP-43
    percent_cells_major_cyt_GFP=num2cell(sum(lookup.MajorCytGFP)/height(lookup)*100); %percentage of cells with majority cytoplasmic TDP-43
    totalp62=height(lookup)-sum(isnan(lookup.Total_p62_puncta)); %total number of cells (without the excluded cells)
    percent_cells_p62_puncta=num2cell((length(find(lookup.Total_p62_puncta))-sum(isnan(lookup.Total_p62_puncta)))/totalp62*100); %percentage of cells with p62 puncta
    T_TGW=[T_TGW;table({str},genotype,percent_cells_major_cyt_TDP,percent_cells_major_cyt_GFP,percent_cells_p62_puncta)];
end
T_TGW.Properties.VariableNames={'mouseID','genotype','percent_cells_MajorCytTDP','percent_cells_MajorCytGFP','percent_cells_p62_puncta'};

for i=1:3
    str=strcat('TGM',num2str(i));
    genotype={'Dync1h1(+/+),TARDBP(-/M)'};
    lookup=TGM(strcmp(TGM.mouseID,str),:);
    percent_cells_major_cyt_TDP=num2cell(sum(lookup.MajorCytTDP)/height(lookup)*100); %percentage of cells with majority cytoplasmic TDP-43
    percent_cells_major_cyt_GFP=num2cell(sum(lookup.MajorCytGFP)/height(lookup)*100); %percentage of cells with majority cytoplasmic TDP-43
    totalp62=height(lookup)-sum(isnan(lookup.Total_p62_puncta)); %total number of cells (without the excluded cells)
    percent_cells_p62_puncta=num2cell((length(find(lookup.Total_p62_puncta))-sum(isnan(lookup.Total_p62_puncta)))/totalp62*100); %percentage of cells with p62 puncta
    T_TGM=[T_TGM;table({str},genotype,percent_cells_major_cyt_TDP,percent_cells_major_cyt_GFP,percent_cells_p62_puncta)];
end
T_TGM.Properties.VariableNames={'mouseID','genotype','percent_cells_MajorCytTDP','percent_cells_MajorCytGFP','percent_cells_p62_puncta'};

for i=1:3
    str=strcat('LWT',num2str(i));
    genotype={'Dync1h1(+/L),TARDBP(-/+)'};
    lookup=LWT(strcmp(LWT.mouseID,str),:);
    percent_cells_major_cyt_TDP=num2cell(sum(lookup.MajorCytTDP)/height(lookup)*100); %percentage of cells with majority cytoplasmic TDP-43
    percent_cells_major_cyt_GFP=num2cell(sum(lookup.MajorCytGFP)/height(lookup)*100); %percentage of cells with majority cytoplasmic TDP-43
    totalp62=height(lookup)-sum(isnan(lookup.Total_p62_puncta)); %total number of cells (without the excluded cells)
    percent_cells_p62_puncta=num2cell((length(find(lookup.Total_p62_puncta))-sum(isnan(lookup.Total_p62_puncta)))/totalp62*100); %percentage of cells with p62 puncta
    T_LWT=[T_LWT;table({str},genotype,percent_cells_major_cyt_TDP,percent_cells_major_cyt_GFP,percent_cells_p62_puncta)];
end
T_LWT.Properties.VariableNames={'mouseID','genotype','percent_cells_MajorCytTDP','percent_cells_MajorCytGFP','percent_cells_p62_puncta'};

for i=1:3
    str=strcat('LMT',num2str(i));
    genotype={'Dync1h1(+/L),TARDBP(-/M)'};
    lookup=LMT(strcmp(LMT.mouseID,str),:);
    percent_cells_major_cyt_TDP=num2cell(sum(lookup.MajorCytTDP)/height(lookup)*100); %percentage of cells with majority cytoplasmic TDP-43
    percent_cells_major_cyt_GFP=num2cell(sum(lookup.MajorCytGFP)/height(lookup)*100); %percentage of cells with majority cytoplasmic TDP-43
    totalp62=height(lookup)-sum(isnan(lookup.Total_p62_puncta)); %total number of cells (without the excluded cells)
    percent_cells_p62_puncta=num2cell((length(find(lookup.Total_p62_puncta))-sum(isnan(lookup.Total_p62_puncta)))/totalp62*100); %percentage of cells with p62 puncta
    T_LMT=[T_LMT;table({str},genotype,percent_cells_major_cyt_TDP,percent_cells_major_cyt_GFP,percent_cells_p62_puncta)];
end
T_LMT.Properties.VariableNames={'mouseID','genotype','percent_cells_MajorCytTDP','percent_cells_MajorCytGFP','percent_cells_p62_puncta'};

T3=[T_TGW;T_TGM;T_LWT;T_LMT];
writetable(T3,'masterfile.xlsx','Sheet','Additional data Batch 2','WriteMode','overwritesheet');
end