clear;
clc;
warning off;
addpath(genpath('./'));

%% dataset
ds={'Dermatology'};
% ds={'Dermatology','Wiki_fea','proteinFold','BDGP_fea','scene','LandUse-21','NUS','CCV','stl10_fea','synthetic'};
dsPath = './dataset/';
resPath = './finalresGB-lmd/';
metric = {'ACC','nmi','Purity','Fscore','Precision','Recall','AR','Entropy'};

%% para setting
lambda1 = 2.^[-2:3:13];
lambda2 = 10.^[-1:1:3];
%% start

for dsi = 1:1:length(ds)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         :length(ds)
    % load data & make folder
    dataName = ds{dsi}; 
    disp(dataName);  % print datasetname
    load(strcat(dsPath,dataName)); 
    k = length(unique(Y));
    numofview = length(X);
    numofsample = length(Y);
    d = (1)*k ;
    matpath = strcat(resPath,dataName);
    txtpath = strcat(resPath,strcat(dataName,'.txt'));
    if (~exist(matpath,'file'))
        mkdir(matpath);
        addpath(genpath(matpath));
    end
    dlmwrite(txtpath, strcat('Dataset:',cellstr(dataName), '  Date:',datestr(now)),'-append','delimiter','','newline','pc');

    %% granular-ball splitting
    fea=cell(numofview,1);
    tempanchor=cell(1,numofview);
    anchor=cell(1,numofview);
    tic;
    for v=1:numofview
        fea{v}=mapstd(X{v}',0,1)';
        anchor{v}=getAnchor(fea{v},numofsample);
    end 
    time_split=toc;
    %% method
    id=1;
    for a = 1:length(lambda2)
        for b=1:length(lambda1)
            tic;
            [U,A,Z,iter,obj,alpha,P] = AWAL_MVGC(X,Y,d,anchor,lambda1(b),lambda2(a));
            U = U ./ repmat(sqrt(sum(U .^ 2, 2)), 1, k);

            [res std] = myNMIACCwithmean(U,Y,k);
            timer(id)  = toc;
            fprintf('lambda2:%d\t lambda1:%d\t Res:%12.6f %12.6f %12.6f %12.6f \tTime:%12.6f \n',[lambda2(a) lambda1(b) res(1) res(2) res(3) res(4) timer(id)]);

            resall{id} = res;
            stdall{id} = std;
            objall{id} = obj;

            dlmwrite(txtpath, [lambda1(b) lambda2(a)  res std timer(id)],'-append','delimiter','\t','newline','pc');
            id=id+1;
            save([resPath,'All_',dataName,'.mat'],'resall','stdall','objall','metric');
        end
    end
    [maxRes,index]=max(cell2mat(resall'));
    filename='result.csv';
    newRow = table({dataName}', maxRes(1), maxRes(2), maxRes(3),maxRes(4),time_split,mean(timer), ...
        'VariableNames', {'Dataset', 'ACC', 'NMI', 'Purity', 'Fscore', 'time_split','Time'});

    if isfile(filename)
        existingData = readtable(filename);
        updatedData = [existingData; newRow];  
        writetable(updatedData, filename);     
    else
        writetable(newRow, filename);
    end
    clear resall stdall objall X Y fea U A Z P 
end


