clc;
clear;

% Read data: The dataset is too large to load once, so train dataset is divided into positive and negative data;

fid = fopen('train_negative.txt'); 
snp_negative = [];
i = 1;
while ~feof(fid) 
    line = fgetl(fid);
    
    if mod(i,3)==2   % read snp data
        line = str2num(line);
        snp_negative = [snp_negative;line];
    end
    
    i = i+1;
end

fclose(fid);

fid = fopen('train_positive.txt'); %read train data
snp_positive = [];
i = 1;
while ~feof(fid) 
    line = fgetl(fid);
    
    if mod(i,3)==2   % read snp 
        line = str2num(line);
        snp_positive = [snp_positive;line];
    end
    
    i = i+1;
end

fclose(fid);

% Add label to snp data;

snp_data=[snp_positive;snp_negative];
snp_label=zeros(297726,1); 

    for i=1:26512
        snp_label(i,1)=1;
    end

    for i=26513:297726
        snp_label(i,1)=0;
    end
    snp0 = [snp_data, snp_label];
    m = randperm(297726);
    snp = snp0(m,:);

data=snp(:,1:51);

label=snp(:,52);

% MRMD selection process;

[feature_select,score]=mrmr_mid_d(data,label,51);
save('feature_select.mat', 'feature_select');

[fea,index]=sort(feature_select);
scores=score(:,index);

%plot scores;
%figure
%plot(scores)


 




