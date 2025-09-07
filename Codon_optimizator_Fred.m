%% The aim of this program is to import a DNA sequence and optimize it using an input codon usage table. Finally, codon usage is compared between both sequences and plotted
% 1. Import a DNA coding sequence as a txt format
% 2. Extract GC content and codon usage of the coding region and save in the 1st tab of an excel file
% 3. Plot codon usage of the input DNA sequence and compare with ribosomal protein usage and HGT cluster protein usage
% 4. Generate graphs to visualize codon usage
% 5. Extract amino acid sequence and save in the 2nd tab of the excel file
% 6. Save imported codon usage table from part 2 (used for graphical comparison of codon usages) in the 3rd tab of the excel file
% 7. Generate a new coding sequence with the given codon usage table and save in the 4th tab of the excel file
% 8. Save new DNA sequence in the 4th tab of the excel file
% 9. Extract GC content and codon usage of the new coding sequence
% 10. Export optimized seq to a new tab of the excel file
% 11. Generate graphs to visualize codon usage
% 12. Plot final protein sequence against input sequence to make sure no AA was changed by mistake 

%% 1. Import a DNA coding sequence as a txt format %%
clear all;clc;close all;font=12;x0=10;y0=10;width=1200;height=600;
[file, folder] = uigetfile({'*.*','All Files (*.*)'}, ...
    'Select a File', ...
    'D:\Dropbox\Boulot Fred\Dirk - HGT and tRNAs\2025 Figures and tables\Github Repositories\Codon optimizator\Example\');
filename=fullfile(folder,file);

%% 2. Extract GC content and codon usage of the coding region and save in the 1st tab of an excel file %%
Freq_inputDNA=Fun_Fred_CodonAnalysis(filename); % Launch function to extract codon usage data from DNA sequence

% Export to excel
[tempDir, tempFile] = fileparts(file); 
filenameexport=fullfile(folder,[tempFile,'_datatable','.xlsx']);
writecell(Freq_inputDNA,filenameexport,'Sheet','native gene');

%% 3. Plot codon usage of the input DNA sequence and compare with ribosomal protein usage and HGT cluster protein usage %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IMPORT CODON USAGE TABLE FOR HGT GENES BASED ON DIRK CLUSTERING DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Import codon usage data using a 2-column datafile (column 1 being codon names and column 2 being codon absolute frequency --> see Fun_Fred_CodonAnalysis_bis function for additional information
% Remark: this uses a similar function as above Fun_Fred_CodonAnalysis. The only difference is that it starts from a codon datatable instead of the raw DNA sequence in a text file.
Dirk_HGTgenecluster = {'Ala-GCA',0.0173300000000000;'Ala-GCC',0.0141400000000000;'Ala-GCG',0.0161800000000000;'Ala-GCT',0.0160400000000000;'Arg-AGA',0.00820000000000000;'Arg-AGG',0.00342000000000000;'Arg-CGA',0.00455000000000000;'Arg-CGC',0.00813000000000000;'Arg-CGG',0.00415000000000000;'Arg-CGT',0.0120000000000000;'Asn-AAC',0.0144900000000000;'Asn-AAT',0.0330600000000000;'Asp-GAC',0.0119400000000000;'Asp-GAT',0.0362900000000000;'Cys-TGC',0.00529000000000000;'Cys-TGT',0.00800000000000000;'Gln-CAA',0.0166700000000000;'Gln-CAG',0.0193800000000000;'Glu-GAA',0.0361400000000000;'Glu-GAG',0.0202700000000000;'Gly-GGA',0.0119800000000000;'Gly-GGC',0.0128200000000000;'Gly-GGG',0.00935000000000000;'Gly-GGT',0.0149800000000000;'His-CAC',0.00565000000000000;'His-CAT',0.0165000000000000;'Ile-ATA',0.0179100000000000;'Ile-ATC',0.0156900000000000;'Ile-ATT',0.0364700000000000;'Leu-CTA',0.00813000000000000;'Leu-CTC',0.00714000000000000;'Leu-CTG',0.0218300000000000;'Leu-CTT',0.0182400000000000;'Leu-TTA',0.0259300000000000;'Leu-TTG',0.0124600000000000;'Lys-AAA',0.0436700000000000;'Lys-AAG',0.0146800000000000;'Met-ATG',0.0253200000000000;'Phe-TTC',0.0103100000000000;'Phe-TTT',0.0330900000000000;'Pro-CCA',0.0100000000000000;'Pro-CCC',0.00505000000000000;'Pro-CCG',0.00870000000000000;'Pro-CCT',0.0106000000000000;'Ser-AGC',0.0108700000000000;'Ser-AGT',0.0133800000000000;'Ser-TCA',0.0133300000000000;'Ser-TCC',0.00682000000000000;'Ser-TCG',0.00639000000000000;'Ser-TCT',0.0116300000000000;'END-TAA',0.00178000000000000;'END-TAG',0;'END-TGA',0;'Thr-ACA',0.0116100000000000;'Thr-ACC',0.0113100000000000;'Thr-ACG',0.0112400000000000;'Thr-ACT',0.0101400000000000;'Trp-TGG',0.0114900000000000;'Tyr-TAC',0.00990000000000000;'Tyr-TAT',0.0254000000000000;'Val-GTA',0.0133300000000000;'Val-GTC',0.00990000000000000;'Val-GTG',0.0114900000000000;'Val-GTT',0.0201300000000000};
Dirk_Ribosomalgenecluster = {'Ala-GCA',0.0198000000000000;'Ala-GCC',0.0137900000000000;'Ala-GCG',0.0277800000000000;'Ala-GCT',0.0288500000000000;'Arg-AGA',0;'Arg-AGG',0;'Arg-CGA',0;'Arg-CGC',0.0172400000000000;'Arg-CGG',0;'Arg-CGT',0.0342900000000000;'Asn-AAC',0.0279700000000000;'Asn-AAT',0.00565000000000000;'Asp-GAC',0.0299100000000000;'Asp-GAT',0.0208300000000000;'Cys-TGC',0.00170000000000000;'Cys-TGT',0;'Gln-CAA',0.00595000000000000;'Gln-CAG',0.0275900000000000;'Glu-GAA',0.0508500000000000;'Glu-GAG',0.0161300000000000;'Gly-GGA',0;'Gly-GGC',0.0341900000000000;'Gly-GGG',0;'Gly-GGT',0.0371900000000000;'His-CAC',0.0119100000000000;'His-CAT',0.00213000000000000;'Ile-ATA',0;'Ile-ATC',0.0370400000000000;'Ile-ATT',0.0169500000000000;'Leu-CTA',0;'Leu-CTC',0;'Leu-CTG',0.0559400000000000;'Leu-CTT',0.00231000000000000;'Leu-TTA',0;'Leu-TTG',0;'Lys-AAA',0.0578500000000000;'Lys-AAG',0.0158700000000000;'Met-ATG',0.0259000000000000;'Phe-TTC',0.0217400000000000;'Phe-TTT',0.00904000000000000;'Pro-CCA',0.00211000000000000;'Pro-CCC',0;'Pro-CCG',0.0223900000000000;'Pro-CCT',0.00263000000000000;'Ser-AGC',0.0117600000000000;'Ser-AGT',0;'Ser-TCA',0;'Ser-TCC',0.0127400000000000;'Ser-TCG',0;'Ser-TCT',0.0168400000000000;'END-TAA',0.00559000000000000;'END-TAG',0;'END-TGA',0;'Thr-ACA',0;'Thr-ACC',0.0277100000000000;'Thr-ACG',0.00462000000000000;'Thr-ACT',0.0166700000000000;'Trp-TGG',0.00281000000000000;'Tyr-TAC',0.0134900000000000;'Tyr-TAT',0.00521000000000000;'Val-GTA',0.0161300000000000;'Val-GTC',0.00990000000000000;'Val-GTG',0.0170500000000000;'Val-GTT',0.0347100000000000};
Freq_Dirk_HGTcluster = Fun_Fred_CodonAnalysis_bis(Dirk_HGTgenecluster);
Freq_Dirk_Ribosomcluster = Fun_Fred_CodonAnalysis_bis(Dirk_Ribosomalgenecluster);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUT CODON TABLE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Freq_new_codons = Freq_Dirk_Ribosomcluster;




%% 4. Generate graphs to visualize codon usage  %%
% Generate bar graphs
% Fun_Fred_Codongraphs(Freq_inputDNA); % Function to generate bar plots
% sgtitle('Input gene codon usage'); % Adds a general title above subplots
% Fun_Fred_Codongraphs(Freq_new_codons); % Function to generate bar plots
% sgtitle('Target codon table'); % Adds a general title above subplots

% Generate scatter plot with regression line
figure;set(gcf,'position',[x0,y0,width,height]);
subplot(1,2,1); % Opti sequence vs ribosomal gene cluster codon usage
Fun_Fred_scatter_reg(cell2mat(Freq_inputDNA(2:62,4)),cell2mat(Freq_new_codons(2:62,4)),'Input gene codon usage','Target codon table',font);
title('Initial vs target table codon usage');
subplot(1,2,2); % Opti sequence vs HGT genes cluster codon usage
Fun_Fred_scatter_reg(cell2mat(Freq_inputDNA(2:62,4)),cell2mat(Freq_Dirk_Ribosomcluster(2:62,4)),'Original DNA sequence','HGT gene cluster',font);
sgtitle('Synonymous codon usage - original DNA sequence');

% Generate codon usage comparison on a straight line
lgd={'Input DNA codon usage', 'Target table codon usage'}; % Legend for raw data graph (second graph)
Fun_Fred_Codontable_comparison(Freq_inputDNA,Freq_new_codons,lgd);
sgtitle('Initial vs target table codon usage');

lgd={'Input DNA codon usage', 'Rib-cluster codon usage'}; % Legend for raw data graph (second graph)
Fun_Fred_Codontable_comparison(Freq_inputDNA,Freq_Dirk_Ribosomcluster,lgd);
sgtitle('Codon usage before optim vs rib-cluster codon usage');

%% 5. Extract amino acid sequence and save in the 2nd tab of the excel file %%
FASTAdata = fastaread((filename), 'blockread', [1 Inf]);
datatext=[FASTAdata.Header FASTAdata.Sequence]; % Full sequence of nucleotides
SeqAA_continuous = aminolookup(nt2aa(datatext)); % Get the amino acid sequence as a continuous text
SeqAA = cellstr(reshape(SeqAA_continuous,3,[])'); % Split into single amino acids
writecell(SeqAA,filenameexport,'Sheet','Protein');

%% 6. Save imported codon usage table from part 2 (used for graphical comparison of codon usages) in the 3rd tab of the excel file %%
if length(Freq_new_codons) > 62
    Freq_new_codons(63:64,:)=[];
end
writecell(Freq_new_codons,filenameexport,'Sheet','Dirk codon table');

%% 7. Generate a new coding sequence with the given codon usage table and save in the 4th tab of the excel file
% Generate sequence of amino acids
SeqAA = nt2aa(datatext,'AlternativeStartCodons',false); % Get continuous 1-letter AA code
SeqAA_count = aacount(SeqAA); % Number of each amino acids - organized in a structure with field name and value

% If a STOP codon is there, eliminate it
if strfind(SeqAA,'*') > 0
    SeqAA(length(SeqAA)) = [];
end

% Generate amino acid table with 1-letter code, 3-letter code and counts of each AA in the protein
Code = fieldnames(SeqAA_count); % Get 1-letter code for each AA
Abbreviation = aminolookup(Code); % Get 3-letter code for each AA
Count = struct2array(SeqAA_count)'; % Get total count of each AA in the protein
idx = 1:20; idx = idx'; % Add index numbers 
AA = table(idx,Code,Abbreviation,Count); % Put data together in a table

% Generate codon table with number of each codons to be used according to new codon usage table and amino acid sequence
codons = Freq_new_codons(:,1); % Get codon names
codons(1,:) = []; % Delete column name to keep only codon names
codons = split(codons,'-');
codons_syn_usage = Freq_new_codons(:,4); codons_syn_usage(1,:) = []; % Format synonymous codon usage same as codon cell array
codons = [codons codons_syn_usage]; % Codon usage table

% Determine number of each codons for each amino acid
for i = 1:size(AA,1)
    j = find(strcmp(codons(:,1), AA.Abbreviation{i}));
    codons(j,4) = num2cell(round(cell2mat(codons(j,3))*AA.Count(i)));
    if sum(cell2mat(codons(j,4))) > AA.Count(i) % If inconsistencies between the sum of codons and the number of AA of one type, because of the round function
        idx = find(cell2mat(codons(j,4)) == max(cell2mat(codons(j,4))));
        idx = j(idx);
        idx = max(idx); % In case there are multiple values
        codons(idx,4) = num2cell(max(cell2mat(codons(j,4)))-1);
    end
    if sum(cell2mat(codons(j,4))) < AA.Count(i) % If inconsistencies between the sum of codons and the number of AA of one type, because of the round function
        idx = find(cell2mat(codons(j,4)) == max(cell2mat(codons(j,4))));
        idx = j(idx);
        idx = max(idx); % In case there are multiple values
        codons(idx,4) = num2cell(max(cell2mat(codons(j,4)))+1);
    end    
end

% Generate vectors for each amino acid containing a pool of codons according to probability to use each codon
codon_pool = {};
for i = 1:size(AA,1)
    range = find(strcmp(codons(:,1), AA.Abbreviation{i}));
    str=[];
    for j = min(range):max(range)
        for k = 1:codons{j,4}
            str = [str codons(j,2)]; % Pull together all codons for a given AA
        end
    end
    codon_pool{i}=str; % Each cell of this array corresponds to an AA and contains all available codons
end

% Generate optimized DNA sequence
optim_seq = [];stop_codon = 'TAA';
for i = 1:size(SeqAA,2)
    j = find(strcmp(SeqAA(i),AA.Code)); % Find position in AA table
    idx = randi(length(codon_pool{1,j}));
    new_codon = cell2mat(codon_pool{1,j}(idx));
    codon_pool{1,j}(idx)=[]; % Delete new codon from pool of codon
    optim_seq = [optim_seq new_codon]; % Concatenate new codon with existing sequence
end
optim_seq = [optim_seq stop_codon]; % Add stop codon

%% 8. Save new DNA sequence in the 4th tab of the excel file %%
writematrix(optim_seq,filenameexport,'Sheet','Optiseq');

%% 9. Extract GC content and codon usage of the new coding sequence %%
% Save sequence as text file
filename_optiseq=fullfile(folder,[tempFile,'_Optiseq.txt']);
fileID = fopen(filename_optiseq,'w');
fprintf(fileID,optim_seq);

% Codon usage analysis
Freq_optiDNA=Fun_Fred_CodonAnalysis(filename_optiseq);

%% 10. Export optimized seq to a new tab of the excel file %%
writecell(Freq_optiDNA,filenameexport,'Sheet','Optim DNA');

%% 11. Generate graphs to visualize codon usage  %%
% Generate bar graphs
% Fun_Fred_Codongraphs(Freq_optiDNA); % Function to generate bar plots
% sgtitle('Absolute codon usage of the optimized sequence'); % Adds a general title above subplots

% Generate scatter plot with regression line
figure;set(gcf,'position',[x0,y0,width,height]);
subplot(1,2,1); % Opti sequence vs ribosomal gene cluster codon usage
Fun_Fred_scatter_reg(cell2mat(Freq_optiDNA(2:62,4)),cell2mat(Freq_Dirk_Ribosomcluster(2:62,4)),'Optimized seq','Ribosomal gene cluster',font);
subplot(1,2,2); % Opti sequence vs HGT genes cluster codon usage
Fun_Fred_scatter_reg(cell2mat(Freq_optiDNA(2:62,4)),cell2mat(Freq_Dirk_HGTcluster(2:62,4)),'Optimized seq','HGT gene cluster',font);
sgtitle('Synonymous codon usage correlations');

% Generate codon usage comparison on a straight line
lgd={'Codon usage after optim', 'Ribosome-cluster codon usage'}; % Legend for raw data graph (second graph)
Fun_Fred_Codontable_comparison(Freq_optiDNA,Freq_Dirk_Ribosomcluster,lgd);
sgtitle('Codon usage after optim vs Ribosome-cluster codon usage');

lgd={'Codon usage after optim', 'HGT-cluster codon usage'}; % Legend for raw data graph (second graph)
Fun_Fred_Codontable_comparison(Freq_optiDNA,Freq_Dirk_HGTcluster,lgd);
sgtitle('Codon usage after optim vs HGT-cluster codon usage');

%% 12. Plot final protein sequence against input sequence to make sure no AA was changed by mistake  %%
AAseq_original = nt2aa(datatext); % Get the 1-letter code amino acid sequence of the non-optimized sequence
AAseq_opti = nt2aa(optim_seq); % Same with optimized sequence

if strcmp(AAseq_original,AAseq_opti) == 1
    disp('CORRECT - INPUT PROTEIN SEQUENCE MATCHES THAT OF THE CODON OPTIMIZED DNA')
%     pause
else
    disp('PROBLEM!!!!!!! INPUT PROTEIN SEQUENCE DOES NOT MATCH THAT OF THE CODON OPTIMIZED DNA')
%     pause
end
