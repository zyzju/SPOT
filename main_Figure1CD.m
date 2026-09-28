%% Figure 1C-D public demonstration
clc; clear; close all;

%% Initial parameters
System.N = 512;                    % SLM pixels per side
System.lambda_mm = 1064e-6;        % wavelength: 1064 nm
System.pixelPitch_mm = 15e-3;      % SLM pixel pitch: 15 um
System.focalLength_mm = 9;         % objective focal length parameter
System.relayMagnification = 2;     % SLM-to-pupil magnification

inputWaist_mm = 2.25;              % Gaussian amplitude waist
targetDiameter_um = 15;            % manuscript sample-plane diameter
samplePlaneScale = 0.5;            % sample/design coordinate calibration
superGaussianOrder = 60;           % target amplitude edge order
iterations = 50;                   % fixed update count used in Figure 1D
gsRandomSeed = 20230907;            % seed for the reproducible GS rerun
useArchivedGSReference = true;     % true reproduces the published GS panel

outputFolder = fullfile(fileparts(mfilename('fullpath')),'outputs');
if ~isfolder(outputFolder), mkdir(outputFolder); end

%% Input Gaussian field and 15-um target
inputWaist_px = inputWaist_mm/System.pixelPitch_mm;
M = floor(10*inputWaist_px)+1;
if mod(M,2)==1, M=M+1; end

xSLM_mm = ((1:M)-M/2)*System.pixelPitch_mm;
[Xslm,Yslm] = meshgrid(xSLM_mm,xSLM_mm);
inputAmplitude = exp(-(Xslm.^2+Yslm.^2)/inputWaist_mm^2);

designPitch_mm = System.focalLength_mm*System.lambda_mm/ ...
    (M*System.pixelPitch_mm*System.relayMagnification);
xSample_mm = ((1:M)-M/2-0.5)*designPitch_mm*samplePlaneScale;
[Xsample,Ysample] = meshgrid(xSample_mm,xSample_mm);
targetRadius_mm = targetDiameter_um/2*1e-3;
targetAmplitude = exp(-((Xsample.^2+Ysample.^2)/targetRadius_mm^2).^ ...
    superGaussianOrder);

%% SPOT calculation (protected implementation)
[spotAmplitude,spotPhaseFull,spotPhaseSLM512] = ...
    SPOT_Template(inputAmplitude,targetAmplitude,inputWaist_mm, ...
    System.pixelPitch_mm,iterations,System.N);
spotMetrics = paperMetrics(spotAmplitude,targetAmplitude,System.N);

%% Conventional GS reference
xSLM512_mm = ((1:System.N)-System.N/2)*System.pixelPitch_mm;
[X512,Y512] = meshgrid(xSLM512_mm,xSLM512_mm);
inputAmplitude512 = exp(-(X512.^2+Y512.^2)/inputWaist_mm^2);
samplePitch512_um = System.focalLength_mm*System.lambda_mm/ ...
    (System.N*System.pixelPitch_mm*System.relayMagnification)* ...
    samplePlaneScale*1e3;
xSample512_um = ((1:System.N)-System.N/2)*samplePitch512_um;
[Xs512,Ys512] = meshgrid(xSample512_um*1e-3,xSample512_um*1e-3);
targetAmplitude512 = exp(-((Xs512.^2+Ys512.^2)/targetRadius_mm^2).^ ...
    superGaussianOrder);

if useArchivedGSReference
    archived = load('GS_Cir_publication_reference.mat','target','phase');
    gsAmplitude = archived.target;
    gsPhase = archived.phase;
    gsLabel = 'GS publication reference';
else
    rng(gsRandomSeed,'twister'); %#ok<UNRCH>
    [gsAmplitude,gsPhase] = conventionalGS(inputAmplitude512, ...
        targetAmplitude512,iterations);
    gsLabel = 'GS fixed-seed rerun';
end
gsMetrics = paperMetrics(gsAmplitude,targetAmplitude512,System.N);

%% Save numerical results
Method = ["SPOT/GSSIP";string(gsLabel)];
U_paper_amplitude = [spotMetrics.U_paper_amplitude; ...
    gsMetrics.U_paper_amplitude];
U_intensity = [spotMetrics.U_intensity;gsMetrics.U_intensity];
RMSE_intensity = [spotMetrics.RMSE;gsMetrics.RMSE];
Eta = [spotMetrics.Eta;gsMetrics.Eta];
Metrics = table(Method,U_paper_amplitude,U_intensity,RMSE_intensity,Eta);
writetable(Metrics,fullfile(outputFolder,'Figure1CD_metrics.csv'));

Parameters = struct('System',System,'InputWaist_mm',inputWaist_mm, ...
    'TargetDiameter_um',targetDiameter_um,'SamplePlaneScale',samplePlaneScale, ...
    'SuperGaussianOrder',superGaussianOrder,'Iterations',iterations, ...
    'CalculationGrid',M,'SamplePitch512_um',samplePitch512_um, ...
    'GSRandomSeed',gsRandomSeed,'UseArchivedGSReference', ...
    useArchivedGSReference);
save(fullfile(outputFolder,'Figure1CD_public_demo_results.mat'), ...
    'spotAmplitude','spotPhaseFull','spotPhaseSLM512','gsAmplitude', ...
    'gsPhase','targetAmplitude','targetAmplitude512','Metrics', ...
    'Parameters','-v7.3');

%% Figure 1C-D manuscript layout
Ispot = abs(spotAmplitude).^2; Ispot = Ispot/max(Ispot(:));
Igs = abs(gsAmplitude).^2; Igs = Igs/max(Igs(:));
center = System.N/2-1;
profileGS = Igs(:,center)/max(Igs(:,center));
profileSPOT = Ispot(:,center)/max(Ispot(:,center));
unwrappedGS = manuscriptUnwrappedPhase(gsPhase,4*pi);
unwrappedSPOT = manuscriptUnwrappedPhase(spotPhaseFull,pi);

fig = figure('Color','w','Position',[60 120 1780 620]);
layout = tiledlayout(fig,2,5,'TileSpacing','compact','Padding','compact');

ax = nexttile(layout,1);
imagesc(ax,mod(gsPhase,2*pi),[0 2*pi]); axis(ax,'image','off');
colormap(ax,hsv(256));
text(ax,-42,256,'GS','Rotation',90,'HorizontalAlignment','center', ...
    'FontWeight','bold','FontSize',10,'Clipping','off');

ax = nexttile(layout,2);
surf(ax,unwrappedGS,'EdgeColor','none'); view(ax,42,32);
formatSurfacePanel(ax,unwrappedGS,4*pi,'4\pi');

ax = nexttile(layout,3);
plotIntensityPanel(ax,xSample512_um,Igs,false);

axD = nexttile(layout,4,[2 2]);
plot(axD,xSample512_um,profileGS,'Color',[0 .65 1], ...
    'LineWidth',1.4,'DisplayName','GS'); hold(axD,'on');
plot(axD,xSample512_um,profileSPOT,'Color',[.9 0 0], ...
    'LineWidth',1.6,'DisplayName','SPOT');
xlim(axD,[-20 20]); ylim(axD,[0 1.05]); box(axD,'on');
xlabel(axD,'x (\mu m)'); ylabel(axD,'Norm. Intensity');
legend(axD,'Location','north','Orientation','horizontal','Box','off');
set(axD,'FontName','Arial','FontSize',10,'LineWidth',1,'TickDir','out');

ax = nexttile(layout,6);
imagesc(ax,mod(spotPhaseSLM512,2*pi),[0 2*pi]); axis(ax,'image','off');
colormap(ax,hsv(256));
text(ax,-42,256,'SPOT','Rotation',90,'HorizontalAlignment','center', ...
    'FontWeight','bold','FontSize',10,'Clipping','off');
phaseBar = colorbar(ax,'southoutside');
phaseBar.Ticks = [0 2*pi]; phaseBar.TickLabels = {'0','2\pi'};
phaseBar.FontSize = 8; phaseBar.Box = 'off';

ax = nexttile(layout,7);
surf(ax,unwrappedSPOT,'EdgeColor','none'); view(ax,42,32);
formatSurfacePanel(ax,unwrappedSPOT,pi,'\pi');

ax = nexttile(layout,8);
plotIntensityPanel(ax,xSample512_um,Ispot,true);

annotation(fig,'textbox',[.004 .935 .025 .04],'String','C', ...
    'LineStyle','none','FontWeight','bold','FontSize',12);
annotation(fig,'textbox',[.596 .94 .025 .04],'String','D', ...
    'LineStyle','none','FontWeight','bold','FontSize',12);
annotation(fig,'rectangle',[.025 .955 .565 .035], ...
    'FaceColor',[.92 .92 .92],'EdgeColor','none');
annotation(fig,'textbox',[.17 .952 .30 .04], ...
    'String','High Uniformity and High Efficiency', ...
    'HorizontalAlignment','center','VerticalAlignment','middle', ...
    'LineStyle','none','FontName','Arial','FontSize',11);
print(fig,fullfile(outputFolder,'Figure1CD_public_demo.png'),'-dpng','-r300');
savefig(fig,fullfile(outputFolder,'Figure1CD_public_demo.fig'));

disp(Metrics);
fprintf('Target diameter: %.3f um\n',targetDiameter_um);
fprintf('SPOT: U=%.9f, RMSE=%.9f, eta=%.9f\n', ...
    spotMetrics.U_paper_amplitude,spotMetrics.RMSE,spotMetrics.Eta);

%% Local functions
function metrics = paperMetrics(field,target,N)
amplitude = abs(field); amplitude = amplitude/max(amplitude(:));
intensity = amplitude.^2;
resizedTarget = imresize(target,[N N]);
flatMask = resizedTarget==1;
signal03 = resizedTarget>.3;
signal01 = resizedTarget>.1;

va = amplitude(flatMask);
vi = intensity(flatMask);
metrics.U_paper_amplitude = 1-(max(va)-min(va))/(max(va)+min(va));
metrics.U_intensity = 1-(max(vi)-min(vi))/(max(vi)+min(vi));
intensity = intensity/mean(intensity(signal03));
metrics.RMSE = sqrt(sum((intensity(signal03)-1).^2)/nnz(signal03));
metrics.Eta = sum(intensity(signal01))/sum(intensity(:));
end

function [targetAmplitude,phase] = conventionalGS(inputAmplitude, ...
    desiredAmplitude,iterations) %#ok<DEFNU>
phase = randn(size(inputAmplitude))*2*pi;
field = inputAmplitude.*exp(1i*phase);
for k = 1:iterations
    focalField = fftshift(fft2(field));
    focalField = sqrt(desiredAmplitude).*exp(1i*angle(focalField));
    slmField = ifft2(ifftshift(focalField));
    field = inputAmplitude.*exp(1i*angle(slmField));
end
phase = mod(angle(slmField),2*pi);
targetAmplitude = abs(fftshift(fft2(inputAmplitude.*exp(1i*phase))));
end

function phaseCrop = manuscriptUnwrappedPhase(phase,displayMaximum)
assert(exist('phaseUnwrap','file')==3, ...
    'phaseUnwrap.mexw64 is required for the Figure 1C Unwrap panels.');
center = floor(size(phase,1)/2);
wrappedCrop = double(phase(center-30:center+30,center-30:center+30));
phaseCrop = phaseUnwrap(wrappedCrop);
phaseCrop = phaseCrop-min(phaseCrop(:));
phaseCrop = phaseCrop/max(phaseCrop(:))*displayMaximum;
end

function plotIntensityPanel(ax,x_um,intensity,showColorbar)
imagesc(ax,x_um,x_um,intensity,[0 1]);
axis(ax,'image','xy','off'); xlim(ax,[-20 20]); ylim(ax,[-20 20]);
colormap(ax,hot(256));
hold(ax,'on');
plot(ax,[-20 20],[0 0],'w--','LineWidth',1.1);
if showColorbar
    plot(ax,[-17 -12],[-17 -17],'w-','LineWidth',3);
    intensityBar = colorbar(ax,'southoutside');
    intensityBar.Ticks = [0 1]; intensityBar.TickLabels = {'0','1'};
    intensityBar.FontSize = 8; intensityBar.Box = 'off';
end
end

function formatSurfacePanel(ax,surfaceData,maximumValue,maximumLabel)
view(ax,-38,30); axis(ax,'tight','vis3d'); box(ax,'on'); grid(ax,'on');
colormap(ax,jet(256)); caxis(ax,[0 maximumValue]);
set(ax,'XTick',linspace(1,size(surfaceData,2),5), ...
    'YTick',linspace(1,size(surfaceData,1),5), ...
    'ZTick',linspace(0,maximumValue,4), ...
    'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[], ...
    'GridColor',[.72 .72 .72],'GridAlpha',.55,'LineWidth',.7);
phaseBar = colorbar(ax,'eastoutside');
phaseBar.Ticks = [0 maximumValue]; phaseBar.TickLabels = {'0',maximumLabel};
phaseBar.FontSize = 8; phaseBar.Box = 'off';
text(ax,8,size(surfaceData,1)*.83,maximumValue*.04,'Unwrap', ...
    'FontName','Arial','FontSize',10,'Rotation',-22,'Color','k');
end
