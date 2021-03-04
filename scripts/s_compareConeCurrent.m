% s_compareConeCurrent

plotAVG = false;

a =load('/Users/winawerlab/matlab/git/toolboxes/JWLOrientedGabor/data/classification/idealobserver/onlyL/ideal_Classify_coneOutputs_contrast0.1000_pa0_eye00_eccen4.50_defocus0.00_noise-none_sf4.00_lms-1.00.00.0.mat');
b =load('/Users/winawerlab/matlab/git/toolboxes/JWLOrientedGabor/data/classification/idealobserver/onlyL_current/current_Classify_coneOutputs_contrast1.0000_pa0_eye00_eccen4.50_defocus0.00_noise-none_sf4.00_lms-1.00.00.0.mat');
c = load('/Users/winawerlab/matlab/git/toolboxes/JWLOrientedGabor/data/classification/defaultnophaseshift/conecurrent/current_Classify_coneOutputs_contrast1.0000_pa0_eye00_eccen4.50_defocus0.00_noise-random_sf4.00_lms-0.60.30.1.mat');
d = load('/Users/winawerlab/matlab/git/toolboxes/JWLOrientedGabor/data/classification/default/conecurrentRV1/current_Classify_coneOutputs_contrast1.0000_pa0_eye11_eccen4.50_defocus0.00_noise-random_sf4.00_lms-0.60.30.1.mat');

if plotAVG
    f = load('/Volumes/server/Projects/PerformanceFieldsIsetBio/data/classification/conecurrent/conedensity/average_4.5deg_nasal_current/current_Classify_coneOutputs_contrast1.000_pa0_eye11_eccen4.50_defocus0.00_noise-random_sf4.00_AVERAGE.mat');
    g = load('/Volumes/server/Projects/PerformanceFieldsIsetBio/data/classification/conecurrent/conedensity/average_4.5deg_temporal_current/current_Classify_coneOutputs_contrast1.000_pa180_eye11_eccen4.50_defocus0.00_noise-random_sf4.00_AVERAGE.mat');
    h = load('/Volumes/server/Projects/PerformanceFieldsIsetBio/data/classification/conecurrent/conedensity/average_4.5deg_inferior_current/current_Classify_coneOutputs_contrast1.000_pa270_eye11_eccen4.50_defocus0.00_noise-random_sf4.00_AVERAGE.mat');
    j = load('/Volumes/server/Projects/PerformanceFieldsIsetBio/data/classification/conecurrent/conedensity/average_4.5deg_superior_current/current_Classify_coneOutputs_contrast1.000_pa90_eye11_eccen4.50_defocus0.00_noise-random_sf4.00_AVERAGE.mat');
else
    f = load('/Volumes/server/Projects/PerformanceFieldsIsetBio/data/classification/conecurrent/conedensity/run1_4.5deg_nasal_current/current_Classify_coneOutputs_contrast1.000_pa0_eye11_eccen4.50_defocus0.00_noise-random_sf4.00.mat');
    g = load('/Volumes/server/Projects/PerformanceFieldsIsetBio/data/classification/conecurrent/conedensity/run1_4.5deg_temporal_current/current_Classify_coneOutputs_contrast1.000_pa0_eye11_eccen4.50_defocus0.00_noise-random_sf4.00.mat');
    h = load('/Volumes/server/Projects/PerformanceFieldsIsetBio/data/classification/conecurrent/conedensity/run1_4.5deg_inferior_current/current_Classify_coneOutputs_contrast1.000_pa0_eye11_eccen4.50_defocus0.00_noise-random_sf4.00.mat');
    j = load('/Volumes/server/Projects/PerformanceFieldsIsetBio/data/classification/conecurrent/conedensity/run1_4.5deg_superior_current/current_Classify_coneOutputs_contrast1.000_pa0_eye11_eccen4.50_defocus0.00_noise-random_sf4.00.mat');
end

f.exParams = loadExpParams('conedensity', false);

colors = lines(8);

figure; set(gcf, 'Position', [150,378,1414,420], 'color', 'w')
plot(a.expParams.contrastLevels, a.accuracy, 'o-', 'color', colors(1,:)); hold on;
plot(b.expParams.contrastLevelsPC, b.accuracy, 'o-', 'color', colors(2,:))
plot(c.expParams.contrastLevelsPC, c.accuracy,'o-', 'color', colors(3,:))
plot(d.expParams.contrastLevelsPC, d.accuracy,'o-', 'color', colors(4,:))
if plotAVG
    plot(f.exParams.contrastLevelsPC, f.P_AVG,'o-', 'color', colors(5,:))
    plot(f.exParams.contrastLevelsPC, g.P_AVG,'o-', 'color', colors(6,:))
    plot(f.exParams.contrastLevelsPC, h.P_AVG,'o-', 'color', colors(7,:))
    plot(f.exParams.contrastLevelsPC, j.P_AVG,'o-', 'color', colors(8,:))
else
    plot(f.exParams.contrastLevelsPC, f.P,'o-', 'color', colors(5,:))
    plot(f.exParams.contrastLevelsPC, g.P,'o-', 'color', colors(6,:))
    plot(f.exParams.contrastLevelsPC, h.P,'o-', 'color', colors(7,:))
    plot(f.exParams.contrastLevelsPC, j.P,'o-', 'color', colors(8,:))
end

legend({'Cone absorptions no noise, L-cone only', ...
    'Cone current no noise, L-cone only',...
    'Cone current w/ photon noise, LMS-cones, no eyemov & stim phase shift', ...
    'Cone current w/ photon noise, LMS-cones, with eyemov & stim phase shift', ...
    'Cone current nasal w/ photon noise, LMS-cones, with eyemov & stim phase shift', ...
    'Cone current temporal w/ photon noise, LMS-cones, with eyemov & stim phase shift', ...
    'Cone current inferior w/ photon noise, LMS-cones, with eyemov & stim phase shift', ...
    'Cone current superior w/ photon noise, LMS-cones, with eyemov & stim phase shift'},'Location', 'BestOutside')
legend boxoff

plot(10.^-5,a.accuracy(1), 'o','color', colors(1,:))
plot(10.^-5,b.accuracy(1), 'o','color', colors(2,:))
plot(10.^-5,c.accuracy(1), 'o','color', colors(3,:))
plot(10.^-5,d.accuracy(1), 'o','color', colors(4,:))

if plotAVG
    plot(10.^-5,f.P_AVG(1), 'o','color', colors(5,:))
    plot(10.^-5,g.P_AVG(1), 'o','color', colors(6,:))
    plot(10.^-5,h.P_AVG(1), 'o','color', colors(7,:))
    plot(10.^-5,j.P_AVG(1), 'o','color', colors(8,:))
else
    plot(10.^-5,f.P(1), 'o','color', colors(5,:))
    plot(10.^-5,g.P(1), 'o','color', colors(6,:))
    plot(10.^-5,h.P(1), 'o','color', colors(7,:))
    plot(10.^-5,j.P(1), 'o','color', colors(8,:))
end

set(gca, 'XScale', 'log')
xlabel('Stimulus contrast (fraction)')
ylabel('Accuracy (%)')
set(gca, 'TickDir', 'out', 'FontSize',15)
box off

if plotAVG
    title('Psychometric functions for ideal and AVG cone current simulations')
else
    title('Psychometric functions for ideal and single run cone current simulations')
end


