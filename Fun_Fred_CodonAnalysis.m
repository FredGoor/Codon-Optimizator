%% The following function:
% (i)   Imports "datatext", a DNA sequence 
% (ii)  Determines codon usage, GC content, number of codons, amino acids usage and synonymous codon usage
% (iii) Organizes data in a cell array "Freq"

function Freq=Fun_Fred_CodonAnalysis(filename)
FASTAdata = fastaread((filename), 'blockread', [1 Inf]);
datatext=[FASTAdata.Header FASTAdata.Sequence]; % Full sequence of nucleotides
Codons = codoncount(datatext);
NTStruct = basecount(datatext);
GC_content = (NTStruct.C + NTStruct.G)/(NTStruct.C + NTStruct.G + NTStruct.A + NTStruct.T)*100;
T = struct2table(Codons);

% Calculate codon usage data
N_aa = sum(T{:,:},2); % Total number of codons
F = T{:,:}./ N_aa; % Absolute frequency of codon usage

% Annotation vector
AAcode = aminolookup(nt2aa(T.Properties.VariableNames,'AlternativeStartCodons',false)); % Determine amino acid name corresponding to each codon (ex: lys)
Cod_names = append(AAcode.','-',T.Properties.VariableNames); % Paste codon to amino acid (ex: lys-AAA)

% Reorganize and sort data
Freq = [Cod_names',num2cell(F)']; % Codon names in first column and absolute codon usage in second column
Freq = sortrows(Freq,1); % Sort alphabetically to put synonymous codons together
Freq([17,18,19], :) = []; % Eliminate ''END'' codons

% Determine synonymous codon usage in a very inelegant way
Ala=sum(cell2mat(Freq(1:4,2)));   Freq(1:4,3)=num2cell(Ala);
Arg=sum(cell2mat(Freq(5:10,2)));  Freq(5:10,3)=num2cell(Arg);
Asn=sum(cell2mat(Freq(11:12,2))); Freq(11:12,3)=num2cell(Asn);
Asp=sum(cell2mat(Freq(13:14,2))); Freq(13:14,3)=num2cell(Asp);
Cys=sum(cell2mat(Freq(15:16,2))); Freq(15:16,3)=num2cell(Cys);
Gln=sum(cell2mat(Freq(17:18,2))); Freq(17:18,3)=num2cell(Gln);
Glu=sum(cell2mat(Freq(19:20,2))); Freq(19:20,3)=num2cell(Glu);
Gly=sum(cell2mat(Freq(21:24,2))); Freq(21:24,3)=num2cell(Gly);
His=sum(cell2mat(Freq(25:26,2))); Freq(25:26,3)=num2cell(His);
Ile=sum(cell2mat(Freq(27:29,2))); Freq(27:29,3)=num2cell(Ile);
Leu=sum(cell2mat(Freq(30:35,2))); Freq(30:35,3)=num2cell(Leu);
Lys=sum(cell2mat(Freq(36:37,2))); Freq(36:37,3)=num2cell(Lys);
Met=sum(cell2mat(Freq(38,2)));    Freq(38,3)=num2cell(Met);
Phe=sum(cell2mat(Freq(39:40,2))); Freq(39:40,3)=num2cell(Phe);
Pro=sum(cell2mat(Freq(41:44,2))); Freq(41:44,3)=num2cell(Pro);
Ser=sum(cell2mat(Freq(45:50,2))); Freq(45:50,3)=num2cell(Ser);
Thr=sum(cell2mat(Freq(51:54,2))); Freq(51:54,3)=num2cell(Thr);
Trp=sum(cell2mat(Freq(55,2)));    Freq(55,3)=num2cell(Trp);
Tyr=sum(cell2mat(Freq(56:57,2))); Freq(56:57,3)=num2cell(Tyr);
Val=sum(cell2mat(Freq(58:61,2))); Freq(58:61,3)=num2cell(Val);

for i=1:61
    if cell2mat(Freq(i,3))==0 % To avoid NAN
        Freq(i,4)=num2cell(0);
    else
        Freq(i,4)=num2cell(cell2mat(Freq(i,2))/cell2mat(Freq(i,3))); % Compute synonymous codon usage
    end
end

% Add codon number, GC content and column titles
Freq(62,1) = {'Codons'};    Freq(62,2) = num2cell(N_aa);       Freq(62,3) = num2cell(0); Freq(62,4) = num2cell(0); % Total number of codons
Freq(63,1) = {'GC-content'};Freq(63,2) = num2cell(GC_content); Freq(63,3) = num2cell(0); Freq(63,4) = num2cell(0); % GC content
Freq = [{'Codons','Absolute codon usage','Absolute AA usage','Synonymous codon usage'};Freq];




end

