clear all
close all

load('simulationOutput.mat');
DataMatrix = table2array(DataRefThisTime);
[Id,ia,ic] = unique(DataMatrix(:,2));
a_counts = accumarray(ic,1);
IdValueCounts = [Id, a_counts];

DataOrderedInterfered = sortrows(DataMatrix,[2 1 12]);
DataOrderedInterfered(:,1) = floor(DataOrderedInterfered(:,1)) - 10;
DataOrderedInterfered_DirectPath = DataOrderedInterfered;

for i = size(DataOrderedInterfered_DirectPath,1):-1:1
    if DataOrderedInterfered(i,23) ~= -1
        DataOrderedInterfered_DirectPath(i,:) = [];
    end
end


DataOrderedInterfered_DirectPath(:,23:25) = [];
DataAbsoluteInfo_Direct = DataOrderedInterfered_DirectPath;

DataOrderedInterfered_DirectPath(:,13) = DataOrderedInterfered_DirectPath(:,13) - DataOrderedInterfered_DirectPath(:,3);
DataOrderedInterfered_DirectPath(:,14) = DataOrderedInterfered_DirectPath(:,14) - DataOrderedInterfered_DirectPath(:,4);
DataOrderedInterfered_DirectPath(:,18) = DataOrderedInterfered_DirectPath(:,18) - DataOrderedInterfered_DirectPath(:,8);
DataOrderedInterfered_DirectPath(:,19) = DataOrderedInterfered_DirectPath(:,19) - DataOrderedInterfered_DirectPath(:,9);

DataOrderedInterfered_DirectPath(:,3) = DataOrderedInterfered_DirectPath(:,3) - DataOrderedInterfered_DirectPath(:,3);
DataOrderedInterfered_DirectPath(:,4) = DataOrderedInterfered_DirectPath(:,4) - DataOrderedInterfered_DirectPath(:,4);
DataOrderedInterfered_DirectPath(:,8) = DataOrderedInterfered_DirectPath(:,8) - DataOrderedInterfered_DirectPath(:,8);
DataOrderedInterfered_DirectPath(:,9) = DataOrderedInterfered_DirectPath(:,9) - DataOrderedInterfered_DirectPath(:,9);

DataOrderedInterfered_DirectPath(:,6) = real(DataOrderedInterfered_DirectPath(:,6) .* exp(1j*DataOrderedInterfered_DirectPath(:,5)));
DataOrderedInterfered_DirectPath(:,16) = real(DataOrderedInterfered_DirectPath(:,16) .* exp(1j*DataOrderedInterfered_DirectPath(:,15)));

DataOrderedInterfered_DirectPath(:,16) = DataOrderedInterfered_DirectPath(:,16) - DataOrderedInterfered_DirectPath(:,6);
DataOrderedInterfered_DirectPath(:,6) = DataOrderedInterfered_DirectPath(:,6) - DataOrderedInterfered_DirectPath(:,6);

for m = Id'
    for l = [0 10 20 30 40 50 60 70 80 90]
        clear chosenInterferedId moment chosenIdMomentData interferingInfo outputFile;
        chosenInterferedId = m;
        moment = l;
        k = 1;
        for i = 1:(size(DataOrderedInterfered_DirectPath,1))
            if DataOrderedInterfered_DirectPath(i,1) == moment && DataOrderedInterfered_DirectPath(i,2) == chosenInterferedId
                chosenIdMomentData(k,:) =  DataOrderedInterfered_DirectPath(i,:);
                k = k + 1;
            end
        end
        if exist('chosenIdMomentData','var') == 1
            interferingInfo = chosenIdMomentData(:,18:19);
            interferingInfo = [interferingInfo chosenIdMomentData(:,16) zeros(size(chosenIdMomentData,1),1)];
            outputFile = ['InterferingInfo/interferingPosition_Id' num2str(chosenInterferedId) '_moment' num2str(moment) '.mat'];
            save(outputFile,"interferingInfo","chosenInterferedId","moment");
        end
    end
end


otherTargetPosition = [];
myMoment = 0;
myId = 336;

while myId ~= DataAbsoluteInfo_Direct(k,2) || myMoment ~= DataAbsoluteInfo_Direct(k,1)
    k = k + 1;
end
x_myId = DataAbsoluteInfo_Direct(k,8);
y_myId = DataAbsoluteInfo_Direct(k,9);

k = 1;
for i = 1:size(DataAbsoluteInfo_Direct)
    if myId ~= DataAbsoluteInfo_Direct(i,2) && myMoment == DataAbsoluteInfo_Direct(i,1)
        otherTargetPosition(k,:) = [DataAbsoluteInfo_Direct(i,8)-x_myId DataAbsoluteInfo_Direct(i,9)-y_myId];
        k = k + 1;
    end
end
save("OtherTargetPosition_336_0.mat","otherTargetPosition")
save("DataAbsoluteInfo_Direct.mat","DataAbsoluteInfo_Direct")








