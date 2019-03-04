%function [ cylStruct ] = ICF( imageFile, mapFile )
%ICF 
%   Code will step through each vein within the sample data vein mask.
%   Each vein of interest should be labelled with a unique number.
%
%   For each vein, cross-sectional slices will be iterated through.
%   For each slice, vein position and radius are estimated
%   Mean radius, slice-radius, and smoothed-slice-radius are used to
%   estimate partial volume correction.
%
%   If progDisplay is true, the process will be shown graphically.
%
%   CylStruct will give the fitted parameters for each vein.


addpath(genpath('Utility'))
addpath(genpath('Sample_Data'))

devPlots = false;
progDisplay = true;

if ~devPlots
% Development inputs
if ~exist('cylMap','var')    
    imageFile = 'Pilot1_QSM.nii.gz';
    cylFile = 'Pilot1_Vein.nii.gz';
    % Development inputs

    image = load_nii(imageFile);
    image = image.img;

    cylMap = load_nii(cylFile);
    cylMap = cylMap.img;
else
    clearvars -except image cylMap devPlots
end

cylROIs = GetROI_ICF(cylMap);

cylStruct(1:numel(cylROIs)) = struct('Id',0,'Orientation',[0 0 0],'ForegroundMu',0,'BackgroundMu',0,'Foreground',[],'Background',[],'Position',[],'Radius',[]);
tic
for i=1:numel(cylROIs)
   if progDisplay
       figure('Units','Normalize','Position',[0 0 1 1])
   end
   
   cylROI = cylROIs(i);
   
   % Cut out region of image
   x = cylROI.Bounds(1,:);
   y = cylROI.Bounds(2,:);
   z = cylROI.Bounds(3,:);
   imageVol = image(x(1):x(2),y(1):y(2),z(1):z(2));
   imageMask = cylMap(x(1):x(2),y(1):y(2),z(1):z(2))==cylROI.Id;
   
   % Organize so that cyl axis is closest to z-axis
   permuteOrder = mod([0 1 2]+cylROI.Dir,3)+1;
   imageVol = permute(imageVol,permuteOrder);
   imageMask = permute(imageMask,permuteOrder);
   
   % Dilate mask in cross-section only
   originalMask = imageMask;
   imageMask = imdilate(imageMask,strel('disk',1));
   
   % Assume perfect alignment with z-axis for first pass of ICF
   [cylPos, ~] = CalcVolGeom_ICF(imageVol,imageMask,[0 0 1], progDisplay);
      
   % Estimate alignment as linear fit to all slices
   px = polyfit(1:size(imageVol,3),cylPos(:,1)',1);
   py = polyfit(1:size(imageVol,3),cylPos(:,2)',1);
   cylZ = [px(1) py(1) 1];
   cylZ = cylZ./sqrt(sum(cylZ.^2));
   
   % Use estimated alignment for second pass
   [cylPos, cylR] = CalcVolGeom_ICF(imageVol,imageMask,cylZ, progDisplay);
      
   % Re-Estimate alignment as linear fit to all slices
   px = polyfit(1:size(imageVol,3),cylPos(:,1)',1);
   py = polyfit(1:size(imageVol,3),cylPos(:,2)',1);
   cylZ = [px(1) py(1) 1];
   cylZ = cylZ./sqrt(sum(cylZ.^2));
      
   % Store output
   cylStruct(i).Id = cylROI.Id;
   cylStruct(i).Orientation(permuteOrder) = cylZ;
   cylStruct(i).Position = [cylPos (0:size(imageVol,3)-1)'];
   
   % Reorder for slice direction
   cylStruct(i).Position(permuteOrder,:) = cylStruct(i).Position([1 2 3],:);
   
   % Correct for location of cylinder in volume
   cylStruct(i).Position(1,:) = cylStruct(i).Position(1,:)+cylROI.Bounds(1,1);
   cylStruct(i).Position(2,:) = cylStruct(i).Position(2,:)+cylROI.Bounds(2,1);
   cylStruct(i).Position(3,:) = cylStruct(i).Position(3,:)+cylROI.Bounds(3,1);
   
   % Store radius estimates per slice, along with smooth and mean
   cylStruct(i).Radius = repmat(mean(cylR),size(cylR));
   cylStruct(i).RawRadius = cylR;
   cylStruct(i).SmoothRadius = smooth(cylR);
   
   % Store maximum intensity of each slice within mask
   cylStruct(i).Max = squeeze(max(max(imageVol.*imageMask,[],2),[],1));
   
   % Analyse results for volume with radius
   cylVol = CalcSlicePV_ICF(size(imageVol), [cylPos (0:size(imageVol,3)-1)'], cylStruct(i).Radius, cylZ);
   p = fit(cylVol(:),double(imageVol(:)),'poly1','Weights',imageMask(:).*0.05+cylVol(:));
   cylStruct(i).ForegroundMu = p.p1+p.p2;
   cylStruct(i).BackgroundMu = p.p2;
   
   % Analyse results for volume with raw radius
   rawCylVol = CalcSlicePV_ICF(size(imageVol), [cylPos (0:size(imageVol,3)-1)'], cylStruct(i).RawRadius, cylZ);
   p = fit(rawCylVol(:),double(imageVol(:)),'poly1','Weights',imageMask(:).*0.05+rawCylVol(:));
   cylStruct(i).RawForegroundMu = p.p1+p.p2;
   cylStruct(i).RawBackgroundMu = p.p2;
   
   % Analyse results for volume with smooth radius
   smoothCylVol = CalcSlicePV_ICF(size(imageVol), [cylPos (0:size(imageVol,3)-1)'], cylStruct(i).SmoothRadius, cylZ);
   p = fit(smoothCylVol(:),double(imageVol(:)),'poly1','Weights',imageMask(:).*0.05+smoothCylVol(:));
   cylStruct(i).SmoothForegroundMu = p.p1+p.p2;
   cylStruct(i).SmoothBackgroundMu = p.p2;
   
   % Analse results per slice
   sliceInd = repmat(1:size(imageVol,3),[size(imageVol,1)*size(imageVol,2) 1]);
   for j=1:size(imageVol,3)
       p = fit(cylVol(sliceInd==j),double(imageVol(sliceInd==j)),'poly1','Weights',imageMask(sliceInd==j).*0.05+cylVol(sliceInd==j));
       cylStruct(i).Foreground(j) = p.p1+p.p2;
       cylStruct(i).Background(j) = p.p2;
       
       p = fit(rawCylVol(sliceInd==j),double(imageVol(sliceInd==j)),'poly1','Weights',imageMask(sliceInd==j).*0.05+rawCylVol(sliceInd==j));
       cylStruct(i).RawForeground(j) = p.p1+p.p2;
       cylStruct(i).RawBackground(j) = p.p2;
       
       p = fit(smoothCylVol(sliceInd==j),double(imageVol(sliceInd==j)),'poly1','Weights',imageMask(sliceInd==j).*0.05+smoothCylVol(sliceInd==j));
       cylStruct(i).SmoothForeground(j) = p.p1+p.p2;
       cylStruct(i).SmoothBackground(j) = p.p2;
   end
   
   if progDisplay
       % Display results
       sliceVol = reshape(imageVol,[size(imageVol,1)*size(imageVol,2) size(imageVol,3)]);
       rawSliceCyl = reshape(rawCylVol,[size(imageVol,1)*size(imageVol,2) size(imageVol,3)]);
       sliceCyl = reshape(cylVol,[size(imageVol,1)*size(imageVol,2) size(imageVol,3)]);
       smoothSliceCyl = reshape(smoothCylVol,[size(imageVol,1)*size(imageVol,2) size(imageVol,3)]);
       sliceMax = max(sliceVol,[],1);
       [~,~,volInd] = ndgrid(1:size(imageVol,1),1:size(imageVol,2),1:size(imageVol,3));       
       for j=1:size(imageVol,3)
           % Plot per slice
           subplot(2,size(imageVol,3),size(imageVol,3)+j)

           % Max
           plot([0 1],[1 1]*cylStruct(i).Max(j))

           hold on
           % Volume vs. Intensity vs. Radius
           p1 = plot(cylVol(cylVol>0 & volInd==j),imageVol(cylVol>0 & volInd==j),'.');
           p2 = plot(rawCylVol(rawCylVol>0 & volInd==j),imageVol(rawCylVol>0 & volInd==j),'x');
           p3 = plot(smoothCylVol(smoothCylVol>0 & volInd==j),imageVol(smoothCylVol>0 & volInd==j),'o');
           
           plot([0 1],[cylStruct(i).Background(j) cylStruct(i).Foreground(j)],'Color',p1.Color)
           plot([0 1],[cylStruct(i).RawBackground(j) cylStruct(i).RawForeground(j)],'Color',p2.Color)
           plot([0 1],[cylStruct(i).SmoothBackground(j) cylStruct(i).SmoothForeground(j)],'Color',p3.Color)

           % Plot corrected intensity
           hold off

           xlim([-.3 1.1])
           legend({'Max','Fg','Fg-Raw','Fg-Smooth'},'Location','SouthEast')
           title(['Radius = ' num2str(cylStruct(i).Radius(j))])
           pause(0.5)
       end
       pause(1)
   end
end
toc
end

if devPlots
% Developer tests
idArray = [];
rArray = [];
icfArray = [];
maxArray = [];
lengthArray = [];

icfMuArray = [];
rMuArray = [];
maxMuArray= [];

icfStdArray = [];
rStdArray = [];
maxStdArray= [];
for i=1:numel(cylStruct)
cyl = cylStruct(i);

idArray(i) = cyl.Id;
rArray(end+(1:numel(cyl.Radius))) = cyl.Radius;
lengthArray(i) = numel(cyl.Radius);
icfArray(end+(1:numel(cyl.Radius))) = cyl.Foreground;
maxArray(end+(1:numel(cyl.Radius))) = cyl.Max;

rMuArray(i) = mean(cyl.Radius);
icfMuArray(i) = cyl.ForegroundMu;
maxMuArray(i) = mean(cyl.Max);

rStdArray(i) = std(cyl.RawRadius);
icfStdArray(i) = std(cyl.Foreground);
maxStdArray(i) = std(cyl.Max);
end
figure('WindowStyle','Docked')
plot(rArray,maxArray-icfArray,'.')
lsline

disp('# | Length | Radius (mu,std) | Max (mu,std) | ICF (mu, std)')
[idArray' lengthArray' rMuArray' rStdArray' icfMuArray' icfStdArray' maxMuArray' maxStdArray']

[r,p] = corr(rArray',maxArray'-icfArray');
disp(['Correlation between Radius and Max-ICF (R,p) = (' num2str(r) ',' num2str(p) ')'])

cnrArray = [];
nArray = [];
r2Array = [];
for i=1:numel(cylROIs)
if progDisplay
figure('Units','Normalize','Position',[0 0 1 1])
end
cylROI = cylROIs(i);
% Cut out region of image
x = cylROI.Bounds(1,:);
y = cylROI.Bounds(2,:);
z = cylROI.Bounds(3,:);
imageVol = image(x(1):x(2),y(1):y(2),z(1):z(2));
imageMask = cylMap(x(1):x(2),y(1):y(2),z(1):z(2))==cylROI.Id;
nArray(i) = std(imageVol(~imageMask));
cnrArray(i) = mean(cylStruct(i).Foreground-cylStruct(i).Background)/nArray(i);
r2Array(i) = cylStruct(i).Radius(1);
end

end



%end

