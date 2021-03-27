clear;clc; close all  

[xeff, yeff, effBoundary] = makePlot(0.29);

effArea = polyarea(xeff,yeff);  

fprintf('The area of the effective airfoil is %s sq ft',effArea);  

% creating .xls file of boundary  

% xlswrite('effBoundary16deg.xls',effBoundary); % creates .xls file of effective boundary coordinates clear;clc  

% [nodenumber   x-coordinate    y-coordinate    u-vel    y-vel] 

function selection(src,event) 

        val = src.Value; 

        str = src.String; 

        disp(val);

        makePlot(val); 

end 

function [xeff, yeff, effBoundary] = makePlot(val) 

    f = figure(1); 

    clf 

    p = uipanel(f,'Position',[0.0 0.0 0.3 0.1]); 

    c = uicontrol(p,'Style','slider'); 

    c.Value = val; 

    c.Callback = @selection; 

    % [nodenumber   x-coordinate    y-coordinate    u-vel    y-vel]  

%     velocityMatrix = dlmread('Velocity.csv'); 

    velocityMatrix = readmatrix('average_flow_field.txt');

    filteredMatrix = velocityMatrix; 

    tic
    for i = 1:1:size(filteredMatrix,1) 

        if velocityMatrix(i,4) < val*30 

           filteredMatrix(i,1) = ( velocityMatrix(i,1) );  

        else  

           filteredMatrix(i,1) = 0;  

        end  

    end
    toc

    c = 1;  
    tic
    while (c <= size(filteredMatrix,1))   

        if filteredMatrix(c,1) == 0  

           rowToDelete = c;  

           filteredMatrix(rowToDelete, :) = [];  

        else    

           c = c + 1;  

        end  

    end
    toc
    
    tic
    DeffectiveCell = mat2cell(filteredMatrix,size(filteredMatrix,1),[1 1 1 1]);  

    scatter(DeffectiveCell{2},DeffectiveCell{3});   % plots filtered coordinates  

    axis equal  

%     xlim([-.5 1.5])  
% 
%     ylim([-.5 .5])  

  

    pbaspect([1 1 1])  

  toc

    % outlining the effective airfoil  

    k = boundary(DeffectiveCell{2},DeffectiveCell{3});  

    hold on; 

    scatter(DeffectiveCell{2}(k),DeffectiveCell{3}(k));  

    % creating data set of boundary  

    xeff = DeffectiveCell{2}(k);  

    yeff = DeffectiveCell{3}(k);  

    effBoundary(:,1) = xeff;  

    effBoundary(:,2) = yeff;  
end 

 