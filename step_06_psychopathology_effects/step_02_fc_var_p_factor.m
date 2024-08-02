clear
clc

root_dir = 'F:/Cui_Lab/Projects/Connectional_Variability_Axis/';
addpath(genpath(root_dir))
working_dir = [root_dir 'step_06_psychopathology_effects/'];
cd(working_dir)

%%
atlas_list = {'schaefer400','glasser360'};
year_list = {'0y','2y'};
term_list = {'general'};
window_list = [400:100:800];
step_list = [50,100];

for atlas_i = 1:length(atlas_list)
    atlas = atlas_list{atlas_i};

    if strcmp(atlas,'schaefer400')
        reg_flag = 1;
    else
        reg_flag = 0;
    end

    for  year_i = 1:length(year_list)
        data_year = year_list{year_i};

        load([root_dir 'step_06_psychopathology_effects/fc/fc_p_factor_' data_year '_' atlas '.mat'])

        for term_i = 1:length(term_list)
            term_label = term_list{term_i};

            for window_i = 1:length(window_list)
                window_length = window_list(window_i);

                for step_i = 1:length(step_list)
                    step_length = step_list(step_i);
                    fc_var_p_factor_sliding_window(fc_p_factor,data_year,window_length,step_length,atlas,term_label,reg_flag)
                end

            end

        end

    end

end

