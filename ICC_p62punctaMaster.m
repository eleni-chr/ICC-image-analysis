function ICC_p62punctaMaster

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
    
    %Load data
    if isfile('merged_p62puncta.xlsx')
        sheets=sheetnames("merged_p62puncta.xlsx");
        if strcmp(sheets(end),'combined')
            data=readtable('merged_p62puncta.xlsx','Sheet','combined');
            
            %Extract identifier info and append to table
            seq=extractAfter(folname,'_');
            ID=extractBefore(seq,'_');
            sampleID=char(zeros(height(data),length(ID)));
            mouseID=char(zeros(height(data),4));
            genotype=char(zeros(height(data),24));
            for i=1:height(data)
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
            data.sampleID=sampleID;
            data.genotype=genotype;
            data.mouseID=mouseID;
            
            image_num=extractAfter(seq,'_');
            if length(image_num)==1
                image_num=strcat('0',image_num);
            end
            im_num=char(zeros(height(data),2));
            for i=1:height(data)
                im_num(i,:)=image_num;
            end
            data.image_num=im_num;

            T=[T;data];
        end
    end
    cd ..
end
cd(parentdir)

if ~isempty(T)
    writetable(T(:,[1:3,6:end]),'p62punctaMaster.xlsx','Sheet','All data','WriteMode','overwritesheet');
    
    %Add additional sheets to the Excel file with separate genotypes
    NTG=T(strcmp(T.genotype,"Dync1h1(+/+),TARDBP(-/-)"),:);
    LOA=T(strcmp(T.genotype,"Dync1h1(+/L),TARDBP(-/-)"),:);

    writetable(NTG,'p62punctaMaster.xlsx','Sheet','NTG','WriteMode','overwritesheet');
    writetable(LOA,'p62punctaMaster.xlsx','Sheet','LOA','WriteMode','overwritesheet');
end

%% Analyse images in Batch 2
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
    
    %Load data
    if isfile('merged_p62puncta.xlsx')
        sheets=sheetnames("merged_p62puncta.xlsx");
        if strcmp(sheets(end),'combined')
            data=readtable('merged_p62puncta.xlsx','Sheet','combined');
            
            %Extract identifier info and append to table
            seq=extractAfter(folname,'_');
            ID=extractBefore(seq,'_');
            sampleID=char(zeros(height(data),length(ID)));
            mouseID=char(zeros(height(data),4));
            genotype=char(zeros(height(data),24));
            for i=1:height(data)
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
            data.sampleID=sampleID;
            data.genotype=genotype;
            data.mouseID=mouseID;
            
            image_num=extractAfter(seq,'_');
            if length(image_num)==1
                image_num=strcat('0',image_num);
            end
            im_num=char(zeros(height(data),2));
            for i=1:height(data)
                im_num(i,:)=image_num;
            end
            data.image_num=im_num;

            T=[T;data];
        end
    end
    cd ..
end
cd(parentdir)

if ~isempty(T)
    writetable(T(:,[1:3,6:end]),'p62punctaMaster.xlsx','Sheet','All data','WriteMode','append');
    
    %Add additional sheets to the Excel file with separate genotypes
    TGW=T(strcmp(T.genotype,"Dync1h1(+/+),TARDBP(-/+)"),:);
    TGM=T(strcmp(T.genotype,"Dync1h1(+/+),TARDBP(-/M)"),:);
    LWT=T(strcmp(T.genotype,"Dync1h1(+/L),TARDBP(-/+)"),:);
    LMT=T(strcmp(T.genotype,"Dync1h1(+/L),TARDBP(-/M)"),:);

    writetable(TGW,'p62punctaMaster.xlsx','Sheet','WT-Tg','WriteMode','overwritesheet');
    writetable(TGM,'p62punctaMaster.xlsx','Sheet','MUT-Tg','WriteMode','overwritesheet');
    writetable(LWT,'p62punctaMaster.xlsx','Sheet','LOA-WT-Tg','WriteMode','overwritesheet');
    writetable(LMT,'p62punctaMaster.xlsx','Sheet','LOA-MUT-Tg','WriteMode','overwritesheet');
end
end