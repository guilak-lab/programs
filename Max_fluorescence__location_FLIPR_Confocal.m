%wtcntl=readmatrix('2020 all lines.xlsx','Sheet','WT control_resp');
wtveh=readmatrix('Analysis for Paper_resp .xlsx','Sheet','WTveh2');
wt205=readmatrix('Analysis for Paper_resp .xlsx','Sheet','WT2052');
%wt1nM=readmatrix('2020 all lines.xlsx','Sheet','WT 1nM101_resp');
%wt1nM205=readmatrix('2020 all lines.xlsx','Sheet','WT 1nM101+205_resp');
wt10nM=readmatrix('Analysis for Paper_resp .xlsx','Sheet','WT1012');
wt10nM205=readmatrix('Analysis for Paper_resp .xlsx','Sheet','WT101+2052');
%wt100nM=readmatrix('2020 all lines.xlsx','Sheet','WT 100nM101_resp');
%wt100nM205=readmatrix('2020 all lines.xlsx','Sheet','WT 100nM101+205_resp');

%vcntl=readmatrix('2020 all lines.xlsx','Sheet','V620I control_resp');
vveh=readmatrix('Analysis for Paper_resp .xlsx','Sheet','V620Iveh');
v205=readmatrix('Analysis for Paper_resp .xlsx','Sheet','V620I205');
%v1nM=readmatrix('2020 all lines.xlsx','Sheet','V620I 1nM101_resp');
%v1nM205=readmatrix('2020 all lines.xlsx','Sheet','V620I 1nM101+205_resp');
v10nM=readmatrix('Analysis for Paper_resp .xlsx','Sheet','V620I101');
v10nM205=readmatrix('Analysis for Paper_resp .xlsx','Sheet','V620I101+205');
%v100nM=readmatrix('2020 all lines.xlsx','Sheet','V620I 100nM101_resp');
%v100nM205=readmatrix('2020 all lines.xlsx','Sheet','V620I 100nM101+205_resp');

%tcntl=readmatrix('2020 all lines.xlsx','Sheet','T89I control_resp');
tveh=readmatrix('Analysis for Paper_resp .xlsx','Sheet','T89Iveh');
t205=readmatrix('Analysis for Paper_resp .xlsx','Sheet','T89I205');
%t1nM=readmatrix('2020 all lines.xlsx','Sheet','T89I 1nM101_resp');
%t1nM205=readmatrix('2020 all lines.xlsx','Sheet','T89I 1nM101+205_resp');
t10nM=readmatrix('Analysis for Paper_resp .xlsx','Sheet','T89I101');
t10nM205=readmatrix('Analysis for Paper_resp .xlsx','Sheet','T89I101+205');
%t100nM=readmatrix('2020 all lines.xlsx','Sheet','T89I 100nM101_resp');
%t100nM205=readmatrix('2020 all lines.xlsx','Sheet','T89I 100nM101+205_resp');

data=cell(4,1);

% data(1)={wtcntl};
data(1)={wtveh};
data(2)={wt205};
% data(3)={wt1nM};
% data(4)={wt1nM205};
data(3)={wt10nM};
data(4)={wt10nM205};
% data(7)={wt100nM};
% data(8)={wt100nM205};

% data(9)={vcntl};
data(5)={vveh};
data(6)={v205};
% data(11)={v1nM};
% data(12)={v1nM205};
data(7)={v10nM};
data(8)={v10nM205};
% data(15)={v100nM};
% data(16)={v100nM205};

% data(19)={tcntl};
data(9)={tveh};
data(10)={t205};
% data(19)={t1nM};
% data(20)={t1nM205};
data(11)={t10nM};
data(12)={t10nM205};
% data(23)={t100nM};
% data(24)={t100nM205};
%%
respondingtime=cell(4,1);

for a=1:4
    mat=data{a};
    resptime=nan(size(mat,2),1);
    responders=sum(mat./300);
    c=1;
    for b=1:size(mat,2)
        if responders(1,b)>=0.25
            resptime(c,1)=find(mat(:,b)==1,1);
            c=c+1;
        else
        end
        respondingtime(a)={resptime};
    end
end

%%

oscillations=cell(1,12);

for a=1:12
    mat=data{a};
    oscil=zeros(size(mat,2),1);
    for b=1:size(mat,2)
        oss=0;
        for c=1:300
            if c==1
                if mat(c,b)==1 && mat(c+1,b)==1
                oss=oss+1;
                end
            elseif c==300
                 if mat(c,b)==1 && mat(c-1,b)==0
                oss=oss+1;
                 end
            else
                if mat(c,b)==1 && mat(c+1,b)==1 && mat(c-1,b)==0
                oss=oss+1;
                end
            end
            oscil(b)=oss; 
        end
        oscillations(a)={oscil};
    end
    
end

%%
% for i=1:400
%     maxlocs(i,1)=max(data(100:400,i));
%     maxlocs(i,2)=find(data(100:400,i)==max(data(100:400,i)),1);
% end

% maxlocs=zeros(400,1);
% responders = sum(data)./300;
% j=1;
% for i=1:400
%     if responders(1,i)>=0.25
%      maxlocs(j,1)=find(data(:,i)==1,1);
%      j=j+1;
%     else
%      maxlocs(j,1)=0;
%      j=j+1;
%     end
% end

