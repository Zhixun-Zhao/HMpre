clc;
clear;

fid = fopen('test_positive.txt');
sequence = [];
site = [];
snp = [];
i = 1;
while ~feof(fid) 
    line = fgetl(fid);
    if mod(i,3)==0   % read sequence
        sequence = [sequence;line];
    end
    
    if mod(i,3)==1   % read site position 
        line = str2num(line);
        site = [site;line];
    end
    
    if mod(i,3)==2   % read snp 
        line = str2num(line);
        snp = [snp;line];
    end
    
    i = i+1;
end

fclose(fid);

str_site=site;
str_sequence=sequence;


% 4-bit Binary Feature
keySet = {'A','C','G','U'}; 
valueSet = {8,4,2,1};
mapObj = containers.Map(keySet,valueSet);
size=size(str_sequence);   
m=size(1);
n=size(2);
D=zeros(m,n);
for i=1:m
    for j=1:n
        D(i,j)=mapObj(str_sequence(i,j));
    end
end
E=dec2bin(D);  
F=cell(m,n);   
for i=1:m
    for j=1:n
        k=size(1)*(j-1)+i; 
        F(i,j)={E(k,:)};
    end
end
dlmwrite('binary_sequence.txt', F)
B=dlmread('binary_sequence.txt');

% transform sequence into opf and add nucleic acid density value 
keySet_0 = {'A','C','G','U'}; 
valueSet_0 = {7,1,4,2};
mapObj_0 = containers.Map(keySet_0,valueSet_0);

D=zeros(m,n);
for i=1:m
    for j=1:n
        D(i,j)=mapObj_0(str_sequence(i,j));
    end
end
E=dec2bin(D); 
F=cell(m,n);  
for i=1:m
    for j=1:n
        k=size(1)*(j-1)+i; 
        F(i,j)={E(k,:)};
    end
end
dlmwrite('opf_sequence.txt', F)
F=dlmread('opf_sequence.txt');
t=4*n;
H=zeros(m,t);
Density=zeros(1,n); % add single nucleic acid density value
Frequency=zeros(m,4);
for i=1:m
    [Density,Frequency(i,:)]=density_function(F(i,:));
    for j=1:51
        k=4*j;
        l=3*j;
        H(i,k-3:k-1)=F(i,l-2:l);
        H(i,k)=Density(j);
    end
end

%2,3-mer feature
keySet_1 = {'AA','AC','AG','AU','CA','CC','CG','CU','GA','GC','GG','GU','UA','UC','UG','UU'};
valueSet_1 = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16];
mapObj_1 = containers.Map(keySet_1,valueSet_1);
 
neighbor_1=zeros(m,n-1);
for i=1:m
    sample_1=str_sequence(i,:);
      for j=1:n-1
      neighbor_1(i,j) = mapObj_1(sample_1(j:j+1));
      end
end

frequency_1=zeros(m,16);
for i=1:m
    for j=1:n-1
        k=neighbor_1(i,j);
        frequency_1(i,k)=frequency_1(i,k)+1;
    end
end
frequency1=frequency_1./(n-1);

keySet_2 = {'AAA','AAC','AAG','AAU','ACA','ACC','ACG','ACU','AGA','AGC','AGG','AGU','AUA','AUC','AUG','AUU','CAA','CAC','CAG','CAU','CCA','CCC','CCG','CCU','CGA','CGC','CGG','CGU','CUA','CUC','CUG','CUU','GAA','GAC','GAG','GAU','GCA','GCC','GCG','GCU','GGA','GGC','GGG','GGU','GUA','GUC','GUG','GUU','UAA','UAC','UAG','UAU','UCA','UCC','UCG','UCU','UGA','UGC','UGG','UGU','UUA','UUC','UUG','UUU',};
valueSet_2 = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64];
mapObj_2 = containers.Map(keySet_2,valueSet_2);

 
neighbor_2=zeros(m,n-2);
for i=1:m
    sample_2=str_sequence(i,:);
      for j=1:n-2
      neighbor_2(i,j) = mapObj_2(sample_2(j:j+2));
      end
end

frequency_2=zeros(m,64);
for i=1:m
    for j=1:n-2
        k=neighbor_2(i,j);
        frequency_2(i,k)=frequency_2(i,k)+1;
    end
end
frequency2=frequency_2./(n-2);

K_mer=[frequency1,frequency2];

% location information(location site, relative location in transcript)
S = str_site;

% flanking window SNP information 
Snp=zeros(m,51);
index=load('snp_index.mat');
in=index.f;
Snp= snp(:,in);


% calculate shannon entropy,En REn IGs
En=zeros(m,7);
for i=1:m
    En(i,1)=Frequency(i,1);
    En(i,2)=Frequency(i,2);
    En(i,3)=Frequency(i,3);
    En(i,4)=Frequency(i,4);
    if En(i,1) ==0
        e11=0;
        e12=0;
    else
        e11=-Frequency(i,1)*log2(Frequency(i,1));
        e12=Frequency(i,1)*log2(Frequency(i,1)*4);
    end
    
    if En(i,2) ==0
        e21=0;
        e22=0;
    else
        e21=-Frequency(i,2)*log2(Frequency(i,2));
        e22=Frequency(i,2)*log2(Frequency(i,2)*4);
    end
    
    if En(i,3) ==0
        e31=0;
        e32=0;
    else
        e31=-Frequency(i,3)*log2(Frequency(i,3));
        e32=Frequency(i,3)*log2(Frequency(i,3)*4);
    end
    
    if En(i,4) ==0
        e41=0;
        e42=0;
    else
        e41=-Frequency(i,4)*log2(Frequency(i,4));
        e42=Frequency(i,4)*log2(Frequency(i,4)*4);
    end
    
    En(i,5)=e11+e21+e31+e41;
    En(i,6)=e12+e22+e32+e42;
    En(i,7)=En(i,5)-En(i,6);
end

% combine sequence and non-sequence features
Feature=[B,H,K_mer,Snp,S,En];  % B:4-bit binary H:CPD K-mer:1,2,3-mer S:transcript location En:Entropy Information

dlmwrite('test_positive_.txt', Feature);

