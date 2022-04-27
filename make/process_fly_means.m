function [] = process_fly_means()
%% SS_open_loop.:
root = 'E:\DATA\Reafferent';
[FILE,PATH] = uigetfile({'*.mat'},'Select data', root, 'MultiSelect', 'on');
FILE = string(FILE);
n_file = length(FILE);

gamma = nan(n_file,1);
ALL = cell(n_file,1);
for n = 1:length(FILE)
    finfo = strsplit(FILE(n), '_');
    gamma(n) = str2double(finfo(3));
    ALL{n} = load(fullfile(PATH,FILE(n)),'FLY','U','N');
end
[gamma, gI] = sort(gamma);
ALL = ALL(gI);

%% Make table
clearvars -except root FILE ALL gamma
clc

n_gamma = length(gamma);

period = ["base", "learn", "relearn"];
n_period = length(period);

metric = ["gain", "phase"];
n_metric = length(metric);

N_fly = cellfun(@(x) x.N.fly, ALL);
n_all = sum(n_gamma*n_period*mean(N_fly));
T = table(nan(n_all,1), categorical(repmat(period(1), [n_all 1])), nan(n_all,15));
T = splitvars(T);
T.Properties.VariableNames = {'gamma','period', 'fly', ...
                            'H_gain', 'H_phase', 'H_compensation_error', ...
                            'Hpred_gain', 'Hpred_phase', 'Hpred_compensation_error', ...
                            'G_gain', 'G_phase', 'H_gain_slope', 'H_phase_slope', ...
                            'H_gain_slope_R2', 'H_gain_slope_P', 'H_phase_slope_R2', 'H_phase_slope_P'};
Data = T;
for m = 1:n_metric
    p = 1;
    for g = 1:n_gamma
        for n = 1:n_period
            for f = 1:N_fly(g)
            	fly_offset = sum(N_fly(1:g-1));
                Data.gamma(p) = gamma(g);
                Data.fly(p) = f + 0*fly_offset;
                Data.period(p) = categorical(period(n));
                
                Data.H_gain(p) = nanmean(ALL{g}.FLY.(period(n)).H(f).gain);
                Data.H_phase(p) = rad2deg(nanmean(ALL{g}.FLY.(period(n)).H(f).phase));
                Data.H_compensation_error(p) = nanmean(ALL{g}.FLY.(period(n)).H(f).compensation_error);
                
                Data.Hpred_gain(p) = nanmean(ALL{g}.FLY.(period(n)).H_prediction(f).gain);
                Data.Hpred_phase(p) = rad2deg(nanmean(ALL{g}.FLY.(period(n)).H_prediction(f).phase));
                Data.Hpred_compensation_error(p) = nanmean(ALL{g}.FLY.(period(n)).H_prediction(f).compensation_error);
                
                Data.G_gain(p) = nanmean(ALL{g}.FLY.(period(n)).G(f).gain);
                Data.G_phase(p) = rad2deg(nanmean(ALL{g}.FLY.(period(n)).G(f).phase));
                
                Data.H_gain_slope(p) = 60*ALL{g}.FLY.(period(n)).H(f).fit_line.gain.slope;
                Data.H_phase_slope(p) = 60*rad2deg(ALL{g}.FLY.(period(n)).H(f).fit_line.phase.slope);
                
                Data.H_gain_slope_R2(p) = ALL{g}.FLY.(period(n)).H(f).fit_line.gain.R2;
                Data.H_phase_slope_R2(p) = ALL{g}.FLY.(period(n)).H(f).fit_line.phase.R2;
                
                Data.H_gain_slope_P(p) = ALL{g}.FLY.(period(n)).H(f).fit_line.gain.P;
                Data.H_phase_slope_P(p) = ALL{g}.FLY.(period(n)).H(f).fit_line.phase.P;
                
                if isnan(Data.H_gain(p))
                    Data.H_gain_slope(p) = nan;
                    Data.H_phase_slope(p) = nan;
                end
                p = p + 1;
            end
        end
    end
end

%% Stats
metric = string(Data.Properties.VariableNames(4:end));
n_metric = length(metric);

P = [];
stats = [];
G = {Data.period, Data.gamma};
for m = 1:n_metric
    [P.(metric(m)),~,stats.(metric(m))] = anovan(Data.(metric(m)), G, 'Display', 'off');
    %[P.(metric(m)),~,stats.(metric(m))] = kruskalwallis(Slope.(metric(m)).slope, G{1});
end
% [results,~,~,gnames] = multcompare(stats.gain, 'Alpha', 0.05, "Dimension",[1 2]);

%% Save
disp('Saving...')
fname = 'combined_gamma';
for n = 1:length(FILE)
    filedata = textscan(FILE{n}, '%s', 'delimiter', '_');
    temp_name = [];
    for k = 3:length(filedata{1})-1
        temp_name = [temp_name '_' char(filedata{1}(k))];
    end
    fname = [fname temp_name];
end
savedir = fullfile(root,'processed');
save(fullfile(savedir, [fname '.mat']), 'Data', 'gamma', 'N_fly', 'P', 'stats', '-v7.3');
disp('Done')

end