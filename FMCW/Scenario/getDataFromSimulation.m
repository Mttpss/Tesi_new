load('simulationOutput.mat');
DataMatrix = table2array(DataRefThisTime);
[Id,ia,ic] = unique(DataMatrix(:,2));
a_counts = accumarray(ic,1);
IdValueCounts = [Id, a_counts];
% DataOrderedInterfered = zeros(size(DataMatrix));

% for i = 1:size(Id)
%     disp(i)
%     for k = 1:size(DataMatrix,1)
%         if DataMatrix(k,2) == Id(i)
%             DataOrderedInterfered(k,:) = DataMatrix(k,:);
%         end
%     end
% end

DataOrderedInterfered = sortrows(DataMatrix,[2 1 12]);
DataOrderedInterfered(:,1) = floor(DataOrderedInterfered(:,1)) - 10;
DataOrderedInterfered_DirectPath = DataOrderedInterfered;

for i = size(DataOrderedInterfered_DirectPath,1):-1:1
    if DataOrderedInterfered(i,23) ~= -1
        DataOrderedInterfered_DirectPath(i,:) = [];
    end
end

DataOrderedInterfered_DirectPath(:,23:25) = [];

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