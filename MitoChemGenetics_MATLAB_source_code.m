%% Determination of Z-scores and deltaZscores from sgRNA read counts
%% Import (normalized) read counts matrix
% The input matrix 'sgRNA_read_counts' can be found in MitoChemGenetics_Table_S2

[numData textData rawData] = xlsread('MitoChemGenetics_Table_S2.xlsx','sgRNA_read_counts');

%% Extract names of sgRNAs

for i = 2:77442;
sgnames{i-1,1} = rawData{i,2};
end
 
%% Convert sgRNA-level data into gene-level data by averaging

genenames = unique(sgnames);

for  i = 1:length(genenames); 
    ind = find(ismember(sgnames,genenames{i}));
    if length(ind)>1;
        geneData(i,:) = mean(numData(ind,:));
    else if length(ind) == 1;
            geneData(i,:) = numData(ind,:);
        end
    end
end

%% Calculate log2 fold changes using gene-level data (Day 15 vs. Day 0)

% DMSO
DMSO(:,1) = geneData(:,8-1)-geneData(:,6-1);
DMSO(:,2) = geneData(:,9-1)-geneData(:,7-1);
m_DMSO =(DMSO(:,1)+DMSO(:,2))/2;

% Met
Met(:,1) = geneData(:,10-1)-geneData(:,5-1);
Met(:,2) = geneData(:,11-1)-geneData(:,7-1);
m_Met =(Met(:,1)+Met(:,2))/2;

% Pier
Pier(:,1) = geneData(:,12-1)-geneData(:,5-1);
Pier(:,2) = geneData(:,13-1)-geneData(:,6-1);
Pier(:,3) = geneData(:,14-1)-geneData(:,7-1);
m_Pier =(Pier(:,1)+Pier(:,2)+Pier(:,3))/3;

% Anti
Anti(:,1) = geneData(:,15-1)-geneData(:,5-1);
Anti(:,2) = geneData(:,16-1)-geneData(:,6-1);
Anti(:,3) = geneData(:,17-1)-geneData(:,7-1);
m_Anti =(Anti(:,1)+Anti(:,2)+Anti(:,3))/3;

% Oligo
Oligo(:,1) = geneData(:,18-1)-geneData(:,5-1);
Oligo(:,2) = geneData(:,19-1)-geneData(:,6-1);
Oligo(:,3) = geneData(:,20-1)-geneData(:,7-1);
m_Oligo =(Oligo(:,1)+Oligo(:,2)+Oligo(:,3))/3;

% AO
AO(:,1) = geneData(:,21-1)-geneData(:,5-1);
AO(:,2) = geneData(:,22-1)-geneData(:,6-1);
AO(:,3) = geneData(:,23-1)-geneData(:,7-1);
m_AO =(AO(:,1)+AO(:,2)+AO(:,3))/3;

% EtBr
EtBr(:,1) = geneData(:,24-1)-geneData(:,5-1);
EtBr(:,2) = geneData(:,25-1)-geneData(:,6-1);
EtBr(:,3) = geneData(:,26-1)-geneData(:,7-1);
m_EtBr =(EtBr(:,1)+EtBr(:,2)+EtBr(:,3))/3;

% Chlor
Chlor(:,1) = geneData(:,27-1)-geneData(:,5-1);
Chlor(:,2) = geneData(:,28-1)-geneData(:,6-1);
Chlor(:,3) = geneData(:,29-1)-geneData(:,7-1);
m_Chlor =(Chlor(:,1)+Chlor(:,2)+Chlor(:,3))/3;

m_LFC = [m_DMSO m_Met m_Pier m_Anti m_Oligo m_AO m_EtBr m_Chlor];

%% Keep only expressed genes for Z-score calculation (log2 FPKM > 0)

ind = find(geneData(:,3)==1);
exp_genenames = genenames(ind);
exp_m_LFC = m_LFC(ind,:);

%% Define Null distribution using non-expressed genes (log2 FPKM = -7)

ind = find(geneData(:,2)==1);
null_d = m_LFC(ind,:);
mean_null = mean(null_d);
std_null = std(null_d);
var_null = var(null_d);

%% Calculate Z-scores as shown in MitoChemGenetics_Table_S1 

for i = 1:8;
    for j = 1:length(exp_m_LFC);
        Z(j,i) = (exp_m_LFC(j,i)-mean_null(i))/std_null(i);
    end
end

% Output matrix:
% The rows are ordered by exp_genenames (11102x1)
% The columns are ordered by [DMSO Met Pier Anti Oligo AO EtBr Chlor]

%% Calculate deltaZscores: (Z_drug - Z_DMSO)/sqrt(2) as shown in MitoChemGenetics_Table_S1 

for i = 1:7
    delta_Z(:,i) = (Z(:,i+1)-Z(:,1))/sqrt(2);
end

% Output matrix:
% The rows are ordered by exp_genenames (11102x1)
% The columns are ordered by [Met Pier Anti Oligo AO EtBr Chlor]

%% Calculate log2 fold changes at sgRNA level (Day 15 vs. Day 0) as shown in MitoChemGenetics_Table_S3 

% DMSO
sgDMSO(:,1) = numData(:,8-1)-numData(:,6-1);
sgDMSO(:,2) = numData(:,9-1)-numData(:,7-1);

% Met
sgMet(:,1) = numData(:,10-1)-numData(:,5-1);
sgMet(:,2) = numData(:,11-1)-numData(:,7-1);

% Pier
sgPier(:,1) = numData(:,12-1)-numData(:,5-1);
sgPier(:,2) = numData(:,13-1)-numData(:,6-1);
sgPier(:,3) = numData(:,14-1)-numData(:,7-1);

% Anti
sgAnti(:,1) = numData(:,15-1)-numData(:,5-1);
sgAnti(:,2) = numData(:,16-1)-numData(:,6-1);
sgAnti(:,3) = numData(:,17-1)-numData(:,7-1);

% Oligo
sgOligo(:,1) = numData(:,18-1)-numData(:,5-1);
sgOligo(:,2) = numData(:,19-1)-numData(:,6-1);
sgOligo(:,3) = numData(:,20-1)-numData(:,7-1);

% AO
sgAO(:,1) = numData(:,21-1)-numData(:,5-1);
sgAO(:,2) = numData(:,22-1)-numData(:,6-1);
sgAO(:,3) = numData(:,23-1)-numData(:,7-1);

% EtBr
sgEtBr(:,1) = numData(:,24-1)-numData(:,5-1);
sgEtBr(:,2) = numData(:,25-1)-numData(:,6-1);
sgEtBr(:,3) = numData(:,26-1)-numData(:,7-1);

% Chlor
sgChlor(:,1) = numData(:,27-1)-numData(:,5-1);
sgChlor(:,2) = numData(:,28-1)-numData(:,6-1);
sgChlor(:,3) = numData(:,29-1)-numData(:,7-1);

sgLFC = [sgDMSO sgMet sgPier sgAnti sgOligo sgAO sgEtBr sgChlor];

% Output matrix:
% The rows are sorted by sgnames (77441x1)
% The columns are sorted by [2xDMSO 2xMet 3xPier 3xAnti 3xOligo 3xAO 3xEtBr 3xChlor]