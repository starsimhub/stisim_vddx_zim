"""
Process IHME data
"""

# %% Imports and settings
import pandas as pd
import pylab as pl
import sciris as sc
import numpy as np


if __name__ == '__main__':

    location = 'zimbabwe'
    load_data = True
    save_totals = True
    plot_data = True
    add_past = False
    add_hiv = False

    load_cs_data = True

    value_vars = [str(i) for i in range(2000, 2051)]

    disease = 'tv'
    disease_mapping = {
        'ct': '395_chlamydial_infection.csv',
        'ng': '396_gonococcal_infection.csv',
        'tv': '397_trichomoniasis.csv',
    }


    # Read in giant IHME file. ***Note, this is not to be committed to the repo!***
    if load_data:
        raw_data = pd.read_csv(f'raw_data/{disease_mapping[disease]}')
        zim_data = raw_data.loc[raw_data.location == location.capitalize()]

        var_mapping = {
            'Deaths': disease+'.new_deaths',
            'Incidence': disease+'.new_infections',
            'Prevalence': disease+'.n_infected',
        }
        dalycols = ['DALYs (Disability-Adjusted Life Years)', 'YLDs (Years Lived with Disability)']
        if disease != 'tv':
            dalycols += ['YLLs (Years of Life Lost)']

        if save_totals:
            zim_data_numbers = zim_data.loc[(zim_data.metric == 'Number')]
            df_melt = zim_data_numbers.melt(id_vars=['sex', 'age', 'measure'], value_vars=value_vars, var_name='year')
            df_both = df_melt.loc[df_melt.sex == 'Both']
            df_group = df_both.groupby(by=['measure', 'year'])['value'].sum().reset_index()
            df_save = df_group.pivot(index='year', columns='measure', values='value')
            df_save = df_save.rename(columns=var_mapping)
            df_save = df_save.drop(columns=dalycols)

            if not add_past:
                df_save.to_csv(f'data/{location}_{disease}_data.csv')

            else:
                # Load HIV data and store population sizes
                hiv_data = pd.read_csv('data/zimbabwe_hiv_data.csv')
                df_save = df_save.reset_index()
                df_save.year = df_save.year.astype(int)
                df_merge = pd.merge(hiv_data, df_save, how='outer', on='year')

                # Load population sizes by age and create proportions
                popsizes = pd.read_csv('data/zimbabwe_popsizes.csv')
                popsizes['prop_0_15'] = popsizes.pop_0_15/popsizes.n_alive
                popsizes['prop_15_49'] = popsizes.pop_15_49/popsizes.n_alive
                popsizes['prop_50_100'] = popsizes.pop_50_100/popsizes.n_alive
                prop_popsizes = popsizes.drop(columns=['Unnamed: 0', 'n_alive', 'pop_0_15', 'pop_15_49', 'pop_50_100'])
                dff = pd.merge(df_merge, prop_popsizes)
                dff['pop_0_15'] = dff['prop_0_15']*dff.n_alive
                dff['pop_15_49'] = dff['prop_15_49']*dff.n_alive
                dff['pop_50_100'] = dff['prop_50_100']*dff.n_alive

                if add_hiv:
                    dff.to_csv(f'data/{location}_data.csv')

        # Extract estimates and save to separate files
        measure_map = {
            'prev': 'Prevalence',
            'deaths': 'Deaths',
            'dalys': 'DALYs (Disability-Adjusted Life Years)',
            'inci': 'Incidence',
        }
        for mkey, mlabel in measure_map.items():
            df = zim_data.loc[(zim_data.measure == mlabel) & (zim_data.metric == 'Number')]
            df.to_csv(f'data/{location}_{disease}_{mkey}.csv')

    # Summarize and plot
    if plot_data:

        to_plot = sc.objdict({'prev': 'burden', 'inci': 'new cases'})
        def group_ages(row):
            if row['age'] < 15:
                return '<15'
            if (row['age'] >= 15) & (row['age'] < 25):
               return '15-24'
            if (row['age'] >= 25) & (row['age'] < 35):
               return '25-34'
            if (row['age'] >= 35) & (row['age'] < 45):
               return '35-44'
            if (row['age'] >= 45) & (row['age'] < 55):
               return '45-54'
            if (row['age'] >= 55):
               return '55+'

        for by_var in ['sex', 'age_group', 'big_age_groups']:

            fig, axes = pl.subplots(1, 2, figsize=(8, 4))
            axes = axes.ravel()

            for pi, pkey, plabel in to_plot.enumitems():

                # Read data and get age/sex groupings
                df_orig = pd.read_csv(f'data/{location}_{disease}_{pkey}.csv')
                df_orig['big_age_groups'] = df_orig.apply(group_ages, axis=1)
                if by_var == 'sex':
                    groups = df_orig.sex.unique()
                elif by_var == 'age_group':
                    groups = df_orig.age_group.unique()
                elif by_var == 'big_age_groups':
                    groups = ['<15', '15-24', '25-34', '35-44', '45-54', '55+']  #df_orig.big_age_groups.unique()

                # Handle colors for age
                if by_var == 'age_group':
                    alabels = dict()
                    for ag in groups:
                        if 'to' in ag: alabels[ag] = int(ag.split(' to ')[0])
                    alabels['<1 year'] = 0
                    alabels['95 plus'] = 95
                    color_list = sc.vectocolor(np.array(list(alabels.values())))
                    colors = {alabel: color_list[ag] for ag, alabel in enumerate(alabels.keys())}
                if by_var == 'big_age_groups':
                    alabels = dict()
                    for ag in groups:
                        if '-' in ag: alabels[ag] = int(ag.split('-')[0])
                    alabels['<15'] = 0
                    alabels['55+'] = 55
                    color_list = sc.vectocolor(np.array(list(alabels.values())))
                    colors = {alabel: color_list[ag] for ag, alabel in enumerate(alabels.keys())}

                elif by_var == 'sex':
                    colors = {'Both': 'k', 'Female': '#ee7989', 'Male': '#4682b4'}

                df_melt = df_orig.melt(id_vars=['sex', 'age_group', 'age', 'big_age_groups'], value_vars=value_vars, var_name='year')
                if by_var in ['age_group', 'big_age_groups']: df_melt = df_melt.loc[df_melt.sex == 'Both']
                df_plot = df_melt.groupby(by=[by_var, 'year'])['value'].sum()

                # Make plot
                ax = axes[pi]
                for group in groups:
                    years = df_plot[group].index.values.astype(int)
                    ax.plot(years, df_plot[group].values, label=group, color=colors[group])

                if by_var in ['sex', 'big_age_groups']: ax.legend(loc='best', frameon=False)
                ax.set_ylabel('Number')
                ax.set_title(f'IHME {disease} estimates: {plabel}')
                ax.set_ylim([0, None])
                sc.SIticks(ax)

            sc.figlayout()
            sc.savefig(f'figures/ihme_{disease}_estimates_{by_var}.png')

