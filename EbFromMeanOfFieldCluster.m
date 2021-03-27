% MATLAB Program to take the average velocity values from a set of exported
% velocity fields then generate an effective body boundary based on this
% averaging.
% Author: Thorne Wolfenbarger
% Date Created: 3/25/21
% Date Updated: 3/25/21

clc; close all;
files = listFilesInDirectory("DataFiles/NACA0015PostStallRe375k");
averaged_flow_field = loadAndAverageFlowFields(files);
writematrix(averaged_flow_field, 'average_flow_field.txt');

plotShrinkFactorOptions(averaged_flow_field);
plotAndSaveContourForFlowField(averaged_flow_field);


function plotAndSaveContourForFlowField(A)
    x = A(:,1);
    y = A(:,2);
    z = A(:,4);

    xv = linspace(min(x),max(x),1001); % ‘x’ Vector For Interpolation
    yv = linspace(min(y),max(y),1001); % ‘y’ Vector For Interpolation
    [X,Y] = ndgrid(xv,yv); % Create Interpolation Grids
    Z = griddata(x, y, z, X, Y); % Interpolated Grid Of ‘tally’ Values

    fig = figure();
    contourf(X,Y,Z, 60,'LineColor','none')
    colorbar('southoutside')
    title("Full Field Contour Plot for Averaged Velocity Field")
    saveas(fig,sprintf("NACA0015_NoLineContour.png"))
    close all

    fig = figure();
    contourf(X,Y,Z, 60,'LineColor','none')
    colorbar('southoutside')
    title("Full Field Contour Plot for Averaged Velocity Field")
    xlim([-0.3 4])
    ylim([-1 1])
    saveas(fig,sprintf("NACA0015_NoLineContourZoomed.png"))
    close all

    fig = figure();
    contourf(X,Y,Z, 60)
    colorbar('southoutside')
    saveas(fig,sprintf("NACA0015_YesLineContour.png"))
    close all
end

function plotShrinkFactorOptions(flow_field)
    minV = 11;
    maxV = 13;

    N_PLOTS = 9;
    indexes = 1:N_PLOTS;
    velocities = (indexes-1)/(N_PLOTS-1);
    for shrinkFactor = linspace(0, 0.5, 5)
        fig = figure('Renderer', 'painters', 'Position',...
                     [200 200 900 1100]);

        for i = indexes
           velocity = velocities(i);
           subplot(9,1,i)
           [X,y] = getBoundary(flow_field, velocity, shrinkFactor);
           X = X(2:length(X));
           y = y(2:length(y));

           hold on
           plot(X,y)

           pbaspect([55 9 1])
           xlim([-0.5, 5])
           ylim([-0.3, 0.6])
           title(sprintf('Velocity = %i m/s', velocity*(maxV-minV)+minV))
        end
        sgtitle(sprintf('EB Geometries with shrink factor of %f', shrinkFactor))
        saveas(fig, sprintf('NACA0015_EB_shrink%0.2f.png', shrinkFactor))
    end
end

function file_list = listFilesInDirectory(base_directory)
    files = dir(base_directory);
    files_struct = struct2cell(files);
    full_file_list = files_struct(1,:);
    
    file_list = [];
    for f = full_file_list
        if (strjoin(f) == ".") || (strjoin(f) == "..")
            continue
        end
       file = strjoin([base_directory '/' f],'');
       file_list = [file_list file];
    end
end

function averagedFlowField = loadAndAverageFlowFields(file_list)
    total_flow_field = 0;
    for file_name = file_list
        if total_flow_field == 0
            total_flow_field = readmatrix(file_name);
        else
            total_flow_field = total_flow_field + readmatrix(file_name);
        end
    end
    
    averagedFlowField = total_flow_field / length(file_list);
end

function [X, y] = getBoundary(velocityMatrix, velocityFilterPercent, shrinkFactor)
A = velocityMatrix;
    minV = 11;
    maxV = 13;

filter_value = minV + velocityFilterPercent*(maxV - minV);
rows_to_keep = A(:,4) < filter_value;
% also discard things that are really far
rows_to_keep = rows_to_keep & (A(:,1) >= -0.1);
rows_to_keep = rows_to_keep & (A(:,2) <= 0.5);
rows_to_keep = rows_to_keep & (A(:,2) >= -0.5);

A = A(rows_to_keep,:);
k = boundary(A(:,1), A(:,2), shrinkFactor);

X=A(k,1);
y=A(k,2);
end