
function  out=Stiffness_Data_Lauren_Program(data)

cd ('/Users/Alireza Savadipour/Desktop/Guilak Rotation')

Directory= ('/Users/Alireza Savadipour/Desktop/Annemarie AFM/1031');
Title=('cell14');
R=2.5/1000000;
v=0.499;
k=10.63;
Date=012121;

% Directory= char(get(data.edit5,'string'));
% Title=char(get(data.edit4,'string'));
% R=str2num(cell2mat(get(data.edit3,'string')));
% v=str2num(cell2mat(get(data.edit2,'string')));
% k=str2num(cell2mat(get(data.edit1,'string')));
% Date=str2num(cell2mat(get(data.edit6,'string')));
% rows=str2num(cell2mat(get(data.edit7,'string')));
% columns=str2num(cell2mat(get(data.edit8,'string')));

Directory=[Directory '/' Title '/'];

Curves = 0; 
RowCount=1;
ColumnCount=0;
CombinedForce=[];
CombinedIndentation=[];
ZeroRow=[];
AverageCFRow=[];
AverageCIRow=[];
AverageCF=[];
AverageCI=[];
CheckForOutliers=1;
ModulusPointsFixed=0;
HeightPointsFixed=0;
Invalids=0;
PercentofMaxIndentationLow=0.25 %#ok<NOPRT>
PercentofMaxIndentationHigh=1 %#ok<NOPRT>
out = struct();

allFiles = dir([Directory,Title, '_data' '/*.ibw']);


NumberofSamples = length(allFiles);

Sample=0;
%rows=sqrt(NumberofSamples);
rows = 20;
%columns=sqrt(NumberofSamples);
columns = 20;
ScanSize=[rows columns];
ElasticMap=zeros(rows,columns);
HeightMap=zeros(rows,columns);
RsqMap=zeros(rows,columns);
List=zeros(NumberofSamples,1);
list_count=1;

while Sample<NumberofSamples % This while can be commented out to see if the code would run
    Sample=Sample+1;
    %waitbar(Sample/NumberofSamples,['Sample ' Sample])
if Sample == 210
    dog=0;
end
    allData=IBWread([Directory Title '_data/' allFiles(Sample).name]);  % Convert IBW files 400*2 matrix (*Raw is the force)
    Defl=allData.y(:,2);             
    Raw=allData.y(:,1);   
     % Define Raw and Defl arrays so that only approach curve is used.
    Raw=Raw+1;
    Defl=Defl+1;
    [a,b]=max(Raw);
    if max(Raw)==1 || max(Defl)==1    % getting rid of all the zeros
        Zeros=find(Raw==1);
        Raw=Raw(1:(Zeros-1));
        Defl=Defl(1:(Zeros-1));
    else
        Raw=Raw(1:b);
        Defl=Defl(1:b);
    end
    
    [z,i]=max(Defl);
    Raw=Raw(1:i)-1;
    Defl=Defl(1:i)-1;
    StartDist = min(Raw);
    RawInitial=Raw-min(Raw);
    DeflInitial=Defl-min(Defl);
    PiezoScanDistance=max(Raw)-min(Raw);
    
    if isempty(Raw)==1 || PiezoScanDistance<1e-6
        E1=0;
        ModelFit=0;
        Force=0;
        IndentationPts=0;
        Invalids=Invalids+1;
        h_total=0;
        d=0;
    else
        
        % Zero Raw and Defl values, omitting the initial 2% of data which often
        % contains oscillating noise.
        NewStart=ceil(.02*length(Raw));
        Raw=Raw(NewStart:length(Raw));
        Defl=Defl(NewStart:length(Defl));
        Raw=Raw-min(Raw);
        Defl=Defl-mean(Defl(1:floor(i/2)));
        %Defl=Defl*CorrectionValue;                      % Correction for InvOLS calibration in fluid
        Delta=Raw-Defl;
        Delta=Delta-min(Delta);
        
        % Flatten approach curve and adjust for slope. Show plot of full approach.
        StartofIndent=find(Defl<.25*max(Defl));
        c=floor(max(StartofIndent)*.5);
        Deltaslope=Delta(1:c);
        Deflslope=Defl(1:floor(c));
        Slope=polyfit(Deltaslope,Deflslope,1); %Fit a polynomial to approach curve
        Deflslopeadj=Defl-Delta*Slope(1); 
        Deflslopeadj=Deflslopeadj-mean(Deflslopeadj(1:c));
        Raw1=Raw;
        
        % Calculate indentation slope of F^2/3 vs. delta
        Delta=Raw1-Deflslopeadj;
        Deltaplot=Delta;
        Force23=(Deflslopeadj*k).^(2/3);   % power of 2/3 is from the Hertz model, "Check darling 2002 paper", sth to do with the sahpe of the indenter so check it with the pyramidal tip
        Force23plot=Force23;
        cutoff=find(Force23>PercentofMaxIndentationLow*max(Force23));
        cutoff1=find(Force23>=PercentofMaxIndentationHigh*max(Force23));
        DeltaFit=Delta(cutoff(1):cutoff1(1));
        Force23Fit=Force23(cutoff(1):cutoff1(1));
        SlopeF=polyfit(DeltaFit,Force23Fit,1);
        
        % Calculate E from slope and back-calculate indentation by determining
        % y-intercept, followed by x-intercept analysis
        E0=(SlopeF(1))^(3/2)*3/4*(1-v^2)/(R^.5);    %Slope with some adjustments (shape matters)
        Delta0=-SlopeF(2)*(3/4*(1-v^2)/(E0*R^.5))^(2/3); % Should this be / not *?... If this is one number after running the code so it is the point that the cantilever hits the sample and if it is not one number so sth else is calculated wrong 
        PiezoMovement=find(Delta>=Delta0);
        Contact=min(PiezoMovement);
        IndentationPts=Delta(PiezoMovement);
        IndentationPts=IndentationPts-min(IndentationPts);
        Force=(Deflslopeadj(PiezoMovement))*k;
        Indented(Sample)=max(IndentationPts);
        Force=Force-min(Force);
        ForceMax(Sample)=max(Force);
        % Calculate distance to surface for first indentation and use to calculate
        % heights for all subsequent curves.
        if isreal(E0)==1
            h_total=RawInitial(length(RawInitial)-length(PiezoMovement));
        else
            h_total=RawInitial(end);
        end
               
        E1 = E0;
        SlopeF1 = SlopeF;
        
        % Thin-layer correction uses max indentation as sole input
        d=max(IndentationPts);
        ThinlayerCorrection=1;%(1-2*a0/pi*(R*d).^0.5/h1+4*a0^2/pi^2*((R*d).^0.5/h1).^2-8/pi^3*(a0^3+4*pi^2/15*b0)*((R*d).^0.5/h1).^3+16*a0/pi^4*(a0^3+3*pi^2/5*b0)*((R*d).^0.5/h1).^4);
        E1=E1/ThinlayerCorrection;
        
        if SlopeF1(1)<0 || length(Contact)<1 || length(RawInitial)==length(PiezoMovement)
            E1=0;
            ModelFit=0;
            Invalids=Invalids+1;
        else
            
            % Make force and model curves using Contact point
            ModelFit=4*E1*R^.5/(3*(1-v.^2))*IndentationPts.^(3/2)*ThinlayerCorrection;       % Includes thin-layer correction that uses max indentation as sole input
            ApproachX=Raw1(1:PiezoMovement(1))-max(Raw1(1:PiezoMovement(1)));
            ApproachY=Deflslopeadj(1:PiezoMovement(1))*k;
            
            if Curves==1
                figure(CombinedIndentationPlots)
                hold on
                plot(ApproachX,ApproachY,'.g')
                plot(IndentationPts, Force, '.')
                plot(IndentationPts, ModelFit, 'r')
                hold off
            end
        end
    end
%    figure
%    plot(Deltaplot,Force23plot)
%    hold on
%    plot(Delta0,Force23(Contact),'ro')
%    legend('Force^{2/3} vs. \Delta','Contact point')
   % Populate elastic moduli map and R^2 map with calculated values
    if RowCount<=rows
        if ColumnCount<=(columns-1)
            ColumnCount=ColumnCount+1;
            ElasticMap(RowCount,ColumnCount)=E1;
            HeightMap(RowCount,ColumnCount)=h_total;
            RsqMap(RowCount,ColumnCount)=1-(norm(Force-ModelFit).^2)/(norm(Force-mean(Force)).^2);
            IndentationMap(RowCount,ColumnCount)=max(IndentationPts);
            IndentMap=d;
            StartDistMap(RowCount,ColumnCount)=StartDist;
        else
            ColumnCount=1;
            RowCount=RowCount+1;
            ElasticMap(RowCount,ColumnCount)=E1;
            HeightMap(RowCount,ColumnCount)=h_total;
            RsqMap(RowCount,ColumnCount)=1-(norm(Force-ModelFit).^2)/(norm(Force-mean(Force)).^2);
            IndentationMap(RowCount,ColumnCount)=max(IndentationPts);
            IndentMap=d;
            StartDistMap(RowCount,ColumnCount)=StartDist;
        end
    end
    List(list_count)=E1;
    list_count=list_count+1;
    
     %   Create matrices containing all Force and Indentation data for plotting a
    %   single curve and determining average values. The average value is based
    %   on averaging -normalized- curves and should not be used as a solution.
    Force1=Force;
    Indentation1=IndentationPts;
    
    if E1~=0
        if Sample==1
            CombinedForce=Force1;
            CombinedIndentation=Indentation1;
            ZeroRow=[ZeroRow -1];
        else
            if size(CombinedForce,1)<size(Force1,1)
                Rows=size(Force1,1)-size(CombinedForce,1);
                for p=1:Rows
                    CombinedForce=[CombinedForce; ZeroRow];
                    CombinedIndentation=[CombinedIndentation; ZeroRow];
                end
            else
                Rows=size(CombinedForce,1)-size(Force1,1);
                for p=1:Rows
                    Force1=[Force1; -1];
                    Indentation1=[Indentation1; -1];
                end
            end
            
            CombinedForce=[CombinedForce Force1];
            CombinedIndentation=[CombinedIndentation Indentation1];
            ZeroRow=[ZeroRow -1];
        end
    end
end


HeightCorrection = StartDistMap - min(min(StartDistMap));
HeightMap = HeightMap;

OutliersRemovedElasticMap=ElasticMap;
OutliersRemovedRsqMap=RsqMap;
out.elasticmap=ElasticMap;
out.rsqd=RsqMap;

% Flatten height map using linear changes to the adjacent corners.
MaxHeight=max(max(HeightMap));
HeightMap=abs(HeightMap-MaxHeight);
SlopeX=((HeightMap(rows,1)-HeightMap(1,1))+(HeightMap(rows,rows)-HeightMap(1,rows)))/rows/2;
SlopeY=((HeightMap(1,columns)-HeightMap(1,1))+(HeightMap(columns,columns)-HeightMap(columns,1)))/columns/2;
FlatHeightMap=zeros(rows,columns);
for i=1:rows
    for j=1:columns
        FlatHeightMap(i,j)=HeightMap(i,j)-i*SlopeX(1)-j*SlopeY(1);
    end
end
HeightMap=FlatHeightMap-min(min(FlatHeightMap));
OutliersRemovedHeightMap=HeightMap;
out.heightmap=HeightMap;
% Save unmodified data
ForceMax=max(ForceMax);
mkdir(Directory);
save([Directory '\' Title '_Stiffness matrix.txt'],'OutliersRemovedElasticMap','-ascii');
save([Directory '\' Title '_Rsq map.txt'],'OutliersRemovedRsqMap','-ascii');
save([Directory '\' Title '_Height matrix.txt'],'OutliersRemovedHeightMap','-ascii');
save([Directory '\' Title '_All moduli list.txt'],'List','-ascii');
save([Directory '\' Title '_Indentation map.txt'],'IndentationMap','-ascii');
save([Directory '\' Title '_ForceMax.txt'],'ForceMax','-ascii');
save([Directory '\' Title '_Deltaplot.txt'],'Deltaplot','-ascii');
save([Directory '\' Title '_Force23plot.txt'],'Force23plot','-ascii');
out=OutliersRemovedElasticMap;
% Calculate average values for force and indentation data with all curves
% starting at their respective zero points. Trailing -1 values are omitted
% when computing average values.
i=1;
for p=1:size(CombinedForce,1)
    for q=1:size(CombinedForce,2)
        if CombinedIndentation(p,q)~=-1
            AverageCFRow(i)=CombinedForce(p,q);
            AverageCIRow(i)=CombinedIndentation(p,q);
            i=i+1;
        end
    end
    i=1;
    AverageCF(p)=mean(AverageCFRow);
    AverageCI(p)=mean(AverageCIRow);
end

% Remove outliers.
ElasticMap1=ElasticMap;
RsqMap1=RsqMap;
HeightMap1=HeightMap;
while CheckForOutliers==1
    
    CheckForOutliers=0;
    
    for i=1:rows
        for j=1:columns
    
% % If the modulus is 2.5 times standard deviation greater or less than the 
% % average surrounding modulus, then that value is removed and replaced by 
% an average of the surrounding cells in the matrix. The same procedure is
% followed for the heights.
            if i==1
                 if j==1
                    AvgMod=mean([ElasticMap(i,j) ElasticMap(i+1,j+1) ElasticMap(i,j+1) ElasticMap(i+1,j)]);
                    StDevMod=std([ElasticMap(i,j) ElasticMap(i+1,j+1) ElasticMap(i,j+1) ElasticMap(i+1,j)]);
                    AvgModH=mean([HeightMap(i,j) HeightMap(i+1,j+1) HeightMap(i,j+1) HeightMap(i+1,j)]);
                    StDevModH=std([HeightMap(i,j) HeightMap(i+1,j+1) HeightMap(i,j+1) HeightMap(i+1,j)]);
                    if ElasticMap(i,j)>(AvgMod+2.5*StDevMod) || ElasticMap(i,j)<(AvgMod-2.5*StDevMod) || ElasticMap(i,j)==1 || isreal(ElasticMap(i,j))==0
                        ElasticMap1(i,j)=(ElasticMap(i+1,j+1)+ElasticMap(i,j+1)+ElasticMap(i+1,j))/3;
                        RsqMap1(i,j)=(RsqMap(i+1,j+1)+RsqMap(i,j+1)+RsqMap(i+1,j))/3; %Tells how well the function fits the data
                        CheckForOutliers=1;
                        ModulusPointsFixed=ModulusPointsFixed+1;
                        OutliersRemovedElasticMap(i,j)=0;
                        OutliersRemovedRsqMap(i,j)=0;
                   elseif HeightMap(i,j)>(AvgModH+2.5*StDevModH) || HeightMap(i,j)<(AvgModH-2.5*StDevModH) || HeightMap(i,j)==0 || isreal(HeightMap(i,j))==0
                        HeightMap1(i,j)=(HeightMap(i+1,j+1)+HeightMap(i,j+1)+HeightMap(i+1,j))/3;
                        CheckForOutliers=1;
                        HeightPointsFixed=HeightPointsFixed+1;
                        OutliersRemovedHeightMap(i,j)=0;
                    end
                 elseif j==columns
                    AvgMod=mean([ElasticMap(i,j) ElasticMap(i,j-1) ElasticMap(i+1,j) ElasticMap(i+1,j-1)]);
                    StDevMod=std([ElasticMap(i,j) ElasticMap(i,j-1) ElasticMap(i+1,j) ElasticMap(i+1,j-1)]);
                    AvgModH=mean([HeightMap(i,j) HeightMap(i,j-1) HeightMap(i+1,j) HeightMap(i+1,j-1)]);
                    StDevModH=std([HeightMap(i,j) HeightMap(i,j-1) HeightMap(i+1,j) HeightMap(i+1,j-1)]);
                    if ElasticMap(i,j)>(AvgMod+2.5*StDevMod) || ElasticMap(i,j)<(AvgMod-2.5*StDevMod) || ElasticMap(i,j)==1 || isreal(ElasticMap(i,j))==0
                        ElasticMap1(i,j)=(ElasticMap(i,j-1)+ElasticMap(i+1,j)+ElasticMap(i+1,j-1))/3;
                        RsqMap1(i,j)=(RsqMap(i,j-1)+RsqMap(i+1,j)+RsqMap(i+1,j-1))/3;  
                        CheckForOutliers=1;
                        ModulusPointsFixed=ModulusPointsFixed+1;
                        OutliersRemovedElasticMap(i,j)=0;
                        OutliersRemovedRsqMap(i,j)=0;
                   elseif HeightMap(i,j)>(AvgModH+2.5*StDevModH) || HeightMap(i,j)<(AvgModH-2.5*StDevModH) || HeightMap(i,j)==0 || isreal(HeightMap(i,j))==0
                        HeightMap1(i,j)=(HeightMap(i,j-1)+HeightMap(i+1,j)+HeightMap(i+1,j-1))/3;
                        CheckForOutliers=1;
                        HeightPointsFixed=HeightPointsFixed+1;
                        OutliersRemovedHeightMap(i,j)=0;
                    end
                 else
                    AvgMod=mean([ElasticMap(i,j) ElasticMap(i,j-1) ElasticMap(i+1,j+1) ElasticMap(i,j+1) ElasticMap(i+1,j) ElasticMap(i+1,j-1)]);
                    StDevMod=std([ElasticMap(i,j) ElasticMap(i,j-1) ElasticMap(i+1,j+1) ElasticMap(i,j+1) ElasticMap(i+1,j) ElasticMap(i+1,j-1)]);
                    AvgModH=mean([HeightMap(i,j) HeightMap(i,j-1) HeightMap(i+1,j+1) HeightMap(i,j+1) HeightMap(i+1,j) HeightMap(i+1,j-1)]);
                    StDevModH=std([HeightMap(i,j) HeightMap(i,j-1) HeightMap(i+1,j+1) HeightMap(i,j+1) HeightMap(i+1,j) HeightMap(i+1,j-1)]);                    
                    if ElasticMap(i,j)>(AvgMod+2.5*StDevMod) || ElasticMap(i,j)<(AvgMod-2.5*StDevMod) || ElasticMap(i,j)==1 || isreal(ElasticMap(i,j))==0
                        ElasticMap1(i,j)=(ElasticMap(i,j-1)+ElasticMap(i+1,j+1)+ElasticMap(i,j+1)+ElasticMap(i+1,j)+ElasticMap(i+1,j-1))/5;
                        RsqMap1(i,j)=(RsqMap(i,j-1)+RsqMap(i+1,j+1)+RsqMap(i,j+1)+RsqMap(i+1,j)+RsqMap(i+1,j-1))/5;
                        CheckForOutliers=1;
                        ModulusPointsFixed=ModulusPointsFixed+1;
                        OutliersRemovedElasticMap(i,j)=0;
                        OutliersRemovedRsqMap(i,j)=0;
                   elseif HeightMap(i,j)>(AvgModH+2.5*StDevModH) || HeightMap(i,j)<(AvgModH-2.5*StDevModH) || HeightMap(i,j)==0 || isreal(HeightMap(i,j))==0
                        HeightMap1(i,j)=(HeightMap(i,j-1)+HeightMap(i+1,j+1)+HeightMap(i,j+1)+HeightMap(i+1,j)+HeightMap(i+1,j-1))/5;
                        CheckForOutliers=1;
                        HeightPointsFixed=HeightPointsFixed+1;
                        OutliersRemovedHeightMap(i,j)=0;
                    end
                 end
            elseif i==rows
                if j==1
                    AvgMod=mean([ElasticMap(i,j) ElasticMap(i-1,j) ElasticMap(i,j+1) ElasticMap(i-1,j+1)]);
                    StDevMod=std([ElasticMap(i,j) ElasticMap(i-1,j) ElasticMap(i,j+1) ElasticMap(i-1,j+1)]);
                    AvgModH=mean([HeightMap(i,j) HeightMap(i-1,j) HeightMap(i,j+1) HeightMap(i-1,j+1)]);
                    StDevModH=std([HeightMap(i,j) HeightMap(i-1,j) HeightMap(i,j+1) HeightMap(i-1,j+1)]);                   
                    if ElasticMap(i,j)>(AvgMod+2.5*StDevMod) || ElasticMap(i,j)<(AvgMod-2.5*StDevMod) || ElasticMap(i,j)==1 || isreal(ElasticMap(i,j))==0
                        ElasticMap1(i,j)=(ElasticMap(i-1,j)+ElasticMap(i,j+1)+ElasticMap(i-1,j+1))/3;
                        RsqMap1(i,j)=(RsqMap(i-1,j)+RsqMap(i,j+1)+RsqMap(i-1,j+1))/3;
                        CheckForOutliers=1;
                        ModulusPointsFixed=ModulusPointsFixed+1;
                        OutliersRemovedElasticMap(i,j)=0;
                        OutliersRemovedRsqMap(i,j)=0;
                   elseif HeightMap(i,j)>(AvgModH+2.5*StDevModH) || HeightMap(i,j)<(AvgModH-2.5*StDevModH) || HeightMap(i,j)==0 || isreal(HeightMap(i,j))==0
                        HeightMap1(i,j)=(HeightMap(i-1,j)+HeightMap(i,j+1)+HeightMap(i-1,j+1))/3;
                        CheckForOutliers=1;
                        HeightPointsFixed=HeightPointsFixed+1;
                        OutliersRemovedHeightMap(i,j)=0;
                    end
                elseif j==columns
                    AvgMod=mean([ElasticMap(i,j) ElasticMap(i-1,j-1) ElasticMap(i,j-1) ElasticMap(i-1,j)]);
                    StDevMod=std([ElasticMap(i,j) ElasticMap(i-1,j-1) ElasticMap(i,j-1) ElasticMap(i-1,j)]);
                    AvgModH=mean([HeightMap(i,j) HeightMap(i-1,j-1) HeightMap(i,j-1) HeightMap(i-1,j)]);
                    StDevModH=std([HeightMap(i,j) HeightMap(i-1,j-1) HeightMap(i,j-1) HeightMap(i-1,j)]);
                    if ElasticMap(i,j)>(AvgMod+2.5*StDevMod) || ElasticMap(i,j)<(AvgMod-2.5*StDevMod) || ElasticMap(i,j)==1 || isreal(ElasticMap(i,j))==0
                        ElasticMap1(i,j)=(ElasticMap(i-1,j-1)+ElasticMap(i,j-1)+ElasticMap(i-1,j))/3;
                        RsqMap1(i,j)=(RsqMap(i-1,j-1)+RsqMap(i,j-1)+RsqMap(i-1,j))/3;
                        CheckForOutliers=1;
                        ModulusPointsFixed=ModulusPointsFixed+1;
                        OutliersRemovedElasticMap(i,j)=0;
                        OutliersRemovedRsqMap(i,j)=0;
                   elseif HeightMap(i,j)>(AvgModH+2.5*StDevModH) || HeightMap(i,j)<(AvgModH-2.5*StDevModH) || HeightMap(i,j)==0 || isreal(HeightMap(i,j))==0
                        HeightMap1(i,j)=(HeightMap(i-1,j-1)+HeightMap(i,j-1)+HeightMap(i-1,j))/3;
                        CheckForOutliers=1;
                        HeightPointsFixed=HeightPointsFixed+1;
                        OutliersRemovedHeightMap(i,j)=0;
                    end
                else
                    AvgMod=mean([ElasticMap(i,j) ElasticMap(i-1,j-1) ElasticMap(i,j-1) ElasticMap(i-1,j) ElasticMap(i,j+1) ElasticMap(i-1,j+1)]);
                    StDevMod=std([ElasticMap(i,j) ElasticMap(i-1,j-1) ElasticMap(i,j-1) ElasticMap(i-1,j) ElasticMap(i,j+1) ElasticMap(i-1,j+1)]);
                    AvgModH=mean([HeightMap(i,j) HeightMap(i-1,j-1) HeightMap(i,j-1) HeightMap(i-1,j) HeightMap(i,j+1) HeightMap(i-1,j+1)]);
                    StDevModH=std([HeightMap(i,j) HeightMap(i-1,j-1) HeightMap(i,j-1) HeightMap(i-1,j) HeightMap(i,j+1) HeightMap(i-1,j+1)]);
                    if ElasticMap(i,j)>(AvgMod+2.5*StDevMod) || ElasticMap(i,j)<(AvgMod-2.5*StDevMod) || ElasticMap(i,j)==1 || isreal(ElasticMap(i,j))==0
                        ElasticMap1(i,j)=(ElasticMap(i-1,j-1)+ElasticMap(i,j-1)+ElasticMap(i-1,j)+ElasticMap(i,j+1)+ElasticMap(i-1,j+1))/5;
                        RsqMap1(i,j)=(RsqMap(i-1,j-1)+RsqMap(i,j-1)+RsqMap(i-1,j)+RsqMap(i,j+1)+RsqMap(i-1,j+1))/5;
                        CheckForOutliers=1;
                        ModulusPointsFixed=ModulusPointsFixed+1;
                        OutliersRemovedElasticMap(i,j)=0;
                        OutliersRemovedRsqMap(i,j)=0;
                   elseif HeightMap(i,j)>(AvgModH+2.5*StDevModH) || HeightMap(i,j)<(AvgModH-2.5*StDevModH) || HeightMap(i,j)==0 || isreal(HeightMap(i,j))==0
                        HeightMap1(i,j)=(HeightMap(i-1,j-1)+HeightMap(i,j-1)+HeightMap(i-1,j)+HeightMap(i,j+1)+HeightMap(i-1,j+1))/5;
                        CheckForOutliers=1;
                        HeightPointsFixed=HeightPointsFixed+1;
                        OutliersRemovedHeightMap(i,j)=0;
                    end
                end
            else
                if j==1
                    AvgMod=mean([ElasticMap(i,j) ElasticMap(i-1,j) ElasticMap(i+1,j+1) ElasticMap(i,j+1) ElasticMap(i+1,j) ElasticMap(i-1,j+1)]);
                    StDevMod=std([ElasticMap(i,j) ElasticMap(i-1,j) ElasticMap(i+1,j+1) ElasticMap(i,j+1) ElasticMap(i+1,j) ElasticMap(i-1,j+1)]);
                    AvgModH=mean([HeightMap(i,j) HeightMap(i-1,j) HeightMap(i+1,j+1) HeightMap(i,j+1) HeightMap(i+1,j) HeightMap(i-1,j+1)]);
                    StDevModH=std([HeightMap(i,j) HeightMap(i-1,j) HeightMap(i+1,j+1) HeightMap(i,j+1) HeightMap(i+1,j) HeightMap(i-1,j+1)]);
                    if ElasticMap(i,j)>(AvgMod+2.5*StDevMod) || ElasticMap(i,j)<(AvgMod-2.5*StDevMod) || ElasticMap(i,j)==1 || isreal(ElasticMap(i,j))==0
                        ElasticMap1(i,j)=(ElasticMap(i-1,j)+ElasticMap(i+1,j+1)+ElasticMap(i,j+1)+ElasticMap(i+1,j)+ElasticMap(i-1,j+1))/5;
                        RsqMap1(i,j)=(RsqMap(i-1,j)+RsqMap(i+1,j+1)+RsqMap(i,j+1)+RsqMap(i+1,j)+RsqMap(i-1,j+1))/5;
                        CheckForOutliers=1;
                        ModulusPointsFixed=ModulusPointsFixed+1;
                        OutliersRemovedElasticMap(i,j)=0;
                        OutliersRemovedRsqMap(i,j)=0;
                   elseif HeightMap(i,j)>(AvgModH+2.5*StDevModH) || HeightMap(i,j)<(AvgModH-2.5*StDevModH) || HeightMap(i,j)==0 || isreal(HeightMap(i,j))==0
                        HeightMap1(i,j)=(HeightMap(i-1,j)+HeightMap(i+1,j+1)+HeightMap(i,j+1)+HeightMap(i+1,j)+HeightMap(i-1,j+1))/5;
                        CheckForOutliers=1;
                        HeightPointsFixed=HeightPointsFixed+1;
                        OutliersRemovedHeightMap(i,j)=0;
                    end
                elseif j==columns
                    AvgMod=mean([ElasticMap(i,j) ElasticMap(i-1,j-1) ElasticMap(i,j-1) ElasticMap(i-1,j) ElasticMap(i+1,j) ElasticMap(i+1,j-1)]);
                    StDevMod=std([ElasticMap(i,j) ElasticMap(i-1,j-1) ElasticMap(i,j-1) ElasticMap(i-1,j) ElasticMap(i+1,j) ElasticMap(i+1,j-1)]);
                    AvgModH=mean([HeightMap(i,j) HeightMap(i-1,j-1) HeightMap(i,j-1) HeightMap(i-1,j) HeightMap(i+1,j) HeightMap(i+1,j-1)]);
                    StDevModH=std([HeightMap(i,j) HeightMap(i-1,j-1) HeightMap(i,j-1) HeightMap(i-1,j) HeightMap(i+1,j) HeightMap(i+1,j-1)]);
                    if ElasticMap(i,j)>(AvgMod+2.5*StDevMod) || ElasticMap(i,j)<(AvgMod-2.5*StDevMod) || ElasticMap(i,j)==1 || isreal(ElasticMap(i,j))==0
                        ElasticMap1(i,j)=(ElasticMap(i-1,j-1)+ElasticMap(i,j-1)+ElasticMap(i-1,j)+ElasticMap(i+1,j)+ElasticMap(i+1,j-1))/5;
                        RsqMap1(i,j)=(RsqMap(i-1,j-1)+RsqMap(i,j-1)+RsqMap(i-1,j)+RsqMap(i+1,j)+RsqMap(i+1,j-1))/5;
                        CheckForOutliers=1;
                        ModulusPointsFixed=ModulusPointsFixed+1;
                        OutliersRemovedElasticMap(i,j)=0;
                        OutliersRemovedRsqMap(i,j)=0;
                   elseif HeightMap(i,j)>(AvgModH+2.5*StDevModH) || HeightMap(i,j)<(AvgModH-2.5*StDevModH) || HeightMap(i,j)==0 || isreal(HeightMap(i,j))==0
                        HeightMap1(i,j)=(HeightMap(i-1,j-1)+HeightMap(i,j-1)+HeightMap(i-1,j)+HeightMap(i+1,j)+HeightMap(i+1,j-1))/5;
                        CheckForOutliers=1;
                        HeightPointsFixed=HeightPointsFixed+1;
                        OutliersRemovedHeightMap(i,j)=0;
                    end
                else
                    AvgMod=mean([ElasticMap(i-1,j-1) ElasticMap(i,j-1) ElasticMap(i-1,j) ElasticMap(i+1,j+1) ElasticMap(i,j+1) ElasticMap(i+1,j) ElasticMap(i+1,j-1) ElasticMap(i-1,j+1)]);
                    StDevMod=std([ElasticMap(i-1,j-1) ElasticMap(i,j-1) ElasticMap(i-1,j) ElasticMap(i+1,j+1) ElasticMap(i,j+1) ElasticMap(i+1,j) ElasticMap(i+1,j-1) ElasticMap(i-1,j+1)]);
                    AvgModH=mean([HeightMap(i-1,j-1) HeightMap(i,j-1) HeightMap(i-1,j) HeightMap(i+1,j+1) HeightMap(i,j+1) HeightMap(i+1,j) HeightMap(i+1,j-1) HeightMap(i-1,j+1)]);
                    StDevModH=std([HeightMap(i-1,j-1) HeightMap(i,j-1) HeightMap(i-1,j) HeightMap(i+1,j+1) HeightMap(i,j+1) HeightMap(i+1,j) HeightMap(i+1,j-1) HeightMap(i-1,j+1)]);
                    if ElasticMap(i,j)>(AvgMod+2.5*StDevMod) || ElasticMap(i,j)<(AvgMod-2.5*StDevMod) || ElasticMap(i,j)==1 || isreal(ElasticMap(i,j))==0
                        ElasticMap1(i,j)=(ElasticMap(i-1,j-1)+ElasticMap(i,j-1)+ElasticMap(i-1,j)+ElasticMap(i+1,j+1)+ElasticMap(i,j+1)+ElasticMap(i+1,j)+ElasticMap(i+1,j-1)+ElasticMap(i-1,j+1))/8;
                        RsqMap1(i,j)=(RsqMap(i-1,j-1)+RsqMap(i,j-1)+RsqMap(i-1,j)+RsqMap(i+1,j+1)+RsqMap(i,j+1)+RsqMap(i+1,j)+RsqMap(i+1,j-1)+RsqMap(i-1,j+1))/8;
                        CheckForOutliers=1;
                        ModulusPointsFixed=ModulusPointsFixed+1;
                        OutliersRemovedElasticMap(i,j)=0;
                        OutliersRemovedRsqMap(i,j)=0;
                   elseif HeightMap(i,j)>(AvgModH+2.5*StDevModH) || HeightMap(i,j)<(AvgModH-2.5*StDevModH) || HeightMap(i,j)==0 || isreal(HeightMap(i,j))==0
                        HeightMap1(i,j)=(HeightMap(i-1,j-1)+HeightMap(i,j-1)+HeightMap(i-1,j)+HeightMap(i+1,j+1)+HeightMap(i,j+1)+HeightMap(i+1,j)+HeightMap(i+1,j-1)+HeightMap(i-1,j+1))/8;
                        CheckForOutliers=1;
                        HeightPointsFixed=HeightPointsFixed+1;
                        OutliersRemovedHeightMap(i,j)=0;
                   end
                end
            end
        end
    end
    ElasticMap=ElasticMap1;
    RsqMap=RsqMap1;
    HeightMap=HeightMap1;
end
Invalids
ModulusPointsFixed
HeightPointsFixed
HeightMap_flat=HeightMap*1e6;

% OutliersRemovedElasticMap

% Create arrays of values associated with outliers removed, cell properties. 
% All non-cell indentations should have E=0, based on the thresholding.
% These arrays will be used to compute average height and modulus values
% for the spread cells.
OutliersRemovedModulusSum=[];
OutliersRemovedRsqSum=[];
OutliersRemovedHeightSum=[];
OutliersRemovedSums=0;
for i=1:rows
    for j=1:columns
        if OutliersRemovedElasticMap(i,j)~=0
            OutliersRemovedSums=OutliersRemovedSums+1;
            OutliersRemovedModulusSum(OutliersRemovedSums)=OutliersRemovedElasticMap(i,j);
            OutliersRemovedRsqSum(OutliersRemovedSums)=OutliersRemovedRsqMap(i,j);
            OutliersRemovedHeightSum(OutliersRemovedSums)=OutliersRemovedHeightMap(i,j);
        end
    end
end

% Calculate average height and modulus, standard deviation, and R^2 of all 
% cells using the original testing data, not the filtered mapping data.
OutliersRemovedAverageElasticModulus=mean2(OutliersRemovedModulusSum)
OutliersRemovedStDev=std2(OutliersRemovedModulusSum);
OutliersRemovedAverageRsq=mean2(OutliersRemovedRsqSum);

% Calculate best fit for the averaged Force vs. Indentation curve. Plot
% averaged data and Hertz model fit.
AverageCF=AverageCF-min(AverageCF);
cd('/Users/Alireza Savadipour/Desktop/Guilak Rotation')

[estimates2, model2] = fitcurveMap(AverageCI,AverageCF,v,R);
model2(estimates2);
Avgfullcurve=4*estimates2(1)*R^(1/2)*AverageCI.^(3/2)/(3*(1-v^2));
 
% figure('Name',[Title ' Indentation Plot'],'NumberTitle','off','Position',[100 100 500 400],'Visible','off')
% IndentationPlot=plot(AverageCI,AverageCF,'.b');
% hold on
% plot(AverageCI,Avgfullcurve,'r');
% ylabel('Force');
% xlabel('Indentation');
% title({sprintf('Elastic indentation (E=%.0f +/- %.0f; R^2=%0.4f)', OutliersRemovedAverageElasticModulus, OutliersRemovedStDev, OutliersRemovedAverageRsq)});%;sprintf('Retraction (E=%.0f +/- %.0f)', RetractionE_average, RetractionE_stdev)});
% hold off

% % Plot all the force-indentation curves on a single figure
% CombinedIndentationPlotExport=figure('Name',[Title ' Combined Indentation Plot'],'NumberTitle','off','Position',[100 100 500 400]);
% 
% CombinedIndentation1=zeros(size(CombinedIndentation,1),size(CombinedIndentation,2));
% CombinedForce1=zeros(size(CombinedForce,1),size(CombinedForce,2));
% for j=1:size(CombinedIndentation,2)
%     for i=1:size(CombinedIndentation,1)
%         if CombinedIndentation(i,j)~=-1 && CombinedForce(i,j)~=-1
%             CombinedIndentation1(i,j)=CombinedIndentation(i,j);
%             CombinedForce1(i,j)=CombinedForce(i,j);
%         end
%     end
% end
% 
% figure(CombinedIndentationPlotExport)
% hold on
% for k=1:size(CombinedIndentation1,1)
%     plot(CombinedIndentation1(k,:),CombinedForce1(k,:),'.b');
% end
% plot(AverageCI,Avgfullcurve,'r');
% ylabel('Force');
% xlabel('Indentation');
% title({sprintf('Elastic indentation (E=%.0f +/- %.0f; R^2=%0.4f)', OutliersRemovedAverageElasticModulus, OutliersRemovedStDev, OutliersRemovedAverageRsq)});%;sprintf('Retraction (E=%.0f +/- %.0f)', RetractionE_average, RetractionE_stdev)});
% hold off

figure('Name',[Title ' Stiffness Map'],'NumberTitle','off','Position',[100 100 500 400]);
[X,Y]=meshgrid(0:ScanSize(1)/(rows-1):ScanSize(2));
AverageHeight2=mean2(OutliersRemovedHeightMap)*1e6;
StDevHeight2=std2(OutliersRemovedHeightMap)*1e6;
SurfaceFlatHeights=surf(X,Y,HeightMap_flat,ElasticMap);
hold on;
xlabel('Width (um)');
ylabel('Length (um)');
zlabel('Height (um)');
hScaleMax=2*max(max(HeightMap_flat));
axis([0 rows 0 columns 0 hScaleMax]);          % Axes scales are set to pre-recorded scansize. Use ZAxisHeight for low profile images.
shading interp
%title(['Leveled Height Map (' num2str(AverageHeight2,'%.2f') ' +/- ' num2str(StDevHeight2,'%.2f') ' um)'; 'Stiffness Map (E=' num2str(OutliersRemovedAverageElasticModulus,'%.2f') ' +/- ' num2str(OutliersRemovedStDev,'%.2f') ', R^2=' num2str(OutliersRemovedAverageRsq,'%1.4f') ', ' num2str(OutliersRemovedSums,'%0.f')]);%; 'ECM = ' num2str(ECM_Modulus/1000,'%.2f') ' kPa (' num2str(ECM_percentage*100,'%0.f') '%), PCM = ' num2str(ECM_Modulus/1000,'%.2f') ' kPa (' num2str(ECM_percentage*100,'%d') '%)']);% 'Stiffness Map (E=%.0f +/- %.0f; R^2=%1.4f; %d points)', OutliersRemovedAverageElasticModulus, OutliersRemovedStDev, OutliersRemovedAverageRsq, OutliersRemovedSums);sprintf('ECM = %.2f kPa (%0.f''%''), PCM = %.2f kPa (%0.f)', ECM_Modulus/1000, ECM_percentage*100, PCM_Modulus/1000, PCM_percentage*100)});      % Include average height measurement in title
title({sprintf('Leveled Height Map (%.2f +/- %.2f um)', AverageHeight2, StDevHeight2);sprintf('Stiffness Map (E=%.0f +/- %.0f; R^2=%1.4f; %d points)', OutliersRemovedAverageElasticModulus, OutliersRemovedStDev, OutliersRemovedAverageRsq, OutliersRemovedSums)});      % Include average height measurement in title
colormap jet
colorbar
hold off
fclose('all')
figure('Name',[Title ' Magnified Stiffness Map'],'NumberTitle','off','Position',[100 100 500 400])
SurfaceFlatHeightsMagnified=surf(X,Y,HeightMap_flat,ElasticMap);
hold on
xlabel('Width (um)');
ylabel('Length (um)');
zlabel('Height (um)');

axis([0 ScanSize(1) 0 ScanSize(2) 0 hScaleMax]);          % Axes scales are set to pre-recorded scansize. Use ZAxisHeight for low profile images.
shading interp
title({sprintf('Magnified, Leveled Height Map (%.2f +/- %.2f um)', AverageHeight2, StDevHeight2);sprintf('Stiffness Map (E=%.0f +/- %.0f; R^2=%1.4f; %d points)', OutliersRemovedAverageElasticModulus, OutliersRemovedStDev, OutliersRemovedAverageRsq, OutliersRemovedSums)});      % Include average height measurement in title
colormap jet
colorbar
hold off

% % Calculate contact radius to look at indentation footprint
% ContactRadius = sqrt(R.*IndentationMap);
% 
% ContactRadiusMap = figure('Name',[Title ' Probe Contact Radius'],'NumberTitle','off','Position',[100 100 500 400]);
% contourf(ContactRadius);
% hold on
% colormap jet
% colorbar
% title('Contact Radius Contour')
% hold off

% % Rsq map figure
% RSQMap = figure('Name',[Title ' R Squared Error'],'NumberTitle','off','Position',[100 100 500 400]);
% contourf(RsqMap);
% hold on
% colormap jet
% colorbar
% title('R Squared Error')
% set(gca,'xtick',[])
% set(gca,'xticklabel',[])
% set(gca,'ytick',[])
% set(gca,'yticklabel',[])
% hold off
% Save the elastic moduli values in 'StiffnessMatrix.dat', R^2 values in 
% 'RsqMap.dat', and the surface plot as 'ElasticMap.jpg'.
saveas(SurfaceFlatHeights,[Directory '\' Title '_Height and stiffness.jpg']);
saveas(SurfaceFlatHeightsMagnified,[Directory '\' Title '_Magnified, height and stiffness.jpg']);
%saveas(IndentationPlot,[Directory '\' Title '_Indentation plot.jpg']);
%saveas(CombinedIndentationPlotExport,[Directory '\' Title '_Combined Indentation plot.jpg']);
%saveas(ContactRadiusMap,[Directory '\' Title '_Contact Radius Map.jpg']);

cd ('/Users/Alireza Savadipour/Desktop/Annemarie AFM/1031/cell14')
Directory1=Directory(1:end);
locations{1}=Directory1;
locations{2}=Title;
locations{3}=Date;

%% Bob edit - 5/14/21
data1.edit5 = Directory;
data1.edit4 = Title;
data1.edit3 = R;
data1.edit2 = v;
data1.edit1 = k;
data1.edit6 = Date;

%[ECM_vec PCM_vec std_error_ECM std_error_PCM]= ECM_PCM(data1);
cd ('C:\Users\Alireza Savadipour\Desktop\Guilak Rotation')
ECM_PCM(data1);

%% revised code

%[ECM_vec PCM_vec std_error_ECM std_error_PCM]= ECM_PCM(out,locations, Directory1,Title);
 save ECM_vec
 save PCM_vec
 save std_error_ECM
 save std_error_PCM
