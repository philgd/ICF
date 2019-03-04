function [ predP, predR ] = CalcVolGeom_ICF( I, mask, predZ, display)
%CALCPVE_VP Summary of this function goes here
%   Detailed explanation goes here

predP = zeros(size(I,3),2);
predR = zeros(size(I,3),1);
bgI = zeros(size(I,3),1);
fgI = zeros(size(I,3),1);

planeDims = size(I);
planeDims(3) = [];

iterMaxN = 15;
stepSize = 1;

for i=1:size(I,3)
        iSlice = reshape(I(:,:,i), planeDims);
        vSlice = reshape(mask(:,:,i), planeDims);
        mSlice = reshape(mask(:,:,i), planeDims);
        
        iSlice = iSlice-min(iSlice(:));

        bgI(i) = mean(iSlice(mSlice==0));
        for j=1:iterMaxN
            % Create Tissue-free Image
            qvMap = mSlice.*max(0,iSlice - stepSize.*bgI(i).*(1-vSlice)); 
            
            % Estimate geometry
            try
                [ predP(i,:), predR(i) ] = CalcSliceGeom_ICF( qvMap, predZ );
            catch ME
                disp(ME)
                disp(ME.message)
                predP(i,:) = [0 0];
                predR(i) = -1;
            end
            
            % Calculate predicted partial volume
            vSlice = CalcVolPV_ICF([planeDims 1],[predP(i,:) 0],predR(i),predZ);
            %[fgI(i), ~, vSlice] = CalcQV_PV( iSlice, predP(i,:), predR(i), predZ);
            
            if display
                % Draw progress
                subplot(2,size(I,3),i)
                imagesc(iSlice,[-0.5 1]*max(abs(I(:))));
                axis image
                colormap gray
                %xlim([1 size(iSlice,2)])
                %ylim([1 size(iSlice,1)])

                % Overlay perimeter
                hold on
                ang = linspace(-pi,pi);
                majorR = predR(i)./predZ(3); % Angle at widest point
                tiltT = atan(predZ(2)./predZ(1)); % Tilt between major axis and image axis
                if isnan(tiltT)
                    tiltT=0;
                end
                rAng = predR(i).*majorR ./ sqrt((predR(i).^2).*(sin(ang-tiltT).^2)+(majorR.^2).*cos(ang-tiltT).^2); % Calculate radius for each angle
                plot(predP(i,2)+0.5, predP(i,1)+0.5,'go')
                plot(predP(i,2)+0.5+(rAng.*cos(ang)),predP(i,1)-(rAng.*sin(ang))+0.5,'g.','LineWidth',3)
                hold off
                drawnow
            end
        end
end

