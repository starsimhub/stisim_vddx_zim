"""
Plot parameter estimates of care seeking alongside prevalence
Requires running run_calibration first to generate the files
'results/zim_sti_calib_{scenario}.obj' for each scenario.
'make_df' creates a big dataframe from the results, and 'plot_pars' plots it.
"""

# Import packages
import sciris as sc
import pylab as pl
import numpy as np
import pandas as pd
from utils import set_font
import seaborn as sns


# %% Run as a script
if __name__ == '__main__':

    to_run = [
        # 'make_df',
        'plot_pars',
    ]

    # Make big dataframe
    scenlabels = {'treat50': 'Treat-half', 'treat80':'Treat-most', 'treat100':'Treat-all'}
    if 'make_df' in to_run:
        dfs = sc.autolist()
        cs_dfs = sc.autolist()  # care seeking for VDS - not by disease
        results = sc.objdict()
        for scenario in ['treat50', 'treat80', 'treat100']:
            calib = sc.loadobj(f'results/zim_sti_calib_{scenario}.obj')
            df = calib.df[:500]
            df['scenario'] = scenlabels[scenario]
            df['ng_p_treat'] = df['ng_p_symp']*df['p_symp_care']*int(scenario.strip('treat'))/100
            df['ct_p_treat'] = df['ct_p_symp']*df['p_symp_care']*int(scenario.strip('treat'))/100
            df['tv_p_treat'] = df['tv_p_symp']*df['p_symp_care']*int(scenario.strip('treat'))/100

            df = df.loc[:, df.columns != 'p_symp_care']
            cs_df = calib.df[:500].loc[:, calib.df.columns.isin(['index', 'p_symp_care'])]
            cs_df['scenario'] = scenario
            dfs += df
            cs_dfs += cs_df

            # Save results
            results[scenario] = calib.sim_results[:50]

        df = pd.concat(dfs)
        cs_df = pd.concat(cs_dfs)

        # Melt dataframe to long form
        vars = ['ng_p_symp', 'ct_p_symp', 'tv_p_symp', 'ng_p_treat', 'ct_p_treat', 'tv_p_treat']
        dfm = df.melt(id_vars=['index', 'mismatch', 'scenario'], value_vars=vars, var_name='variable', value_name='value')

        # Add column for disease
        dfm['disease'] = dfm['variable'].apply(lambda x: x.split('_')[0].upper())
        dfm['par'] = dfm['variable'].apply(lambda x: x[3:])

        # Save
        sc.saveobj('results/zim_sti_calib_df.obj', dfm)
        sc.saveobj('results/zim_sti_care_seeking.obj', cs_df)
        sc.saveobj('results/zim_sti_calib_res.obj', results)

    if 'plot_pars' in to_run: 
        # Load dataframes
        df = sc.loadobj('results/zim_sti_calib_df.obj')
        results = sc.loadobj('results/zim_sti_calib_res.obj')
        # clist = sc.vectocolor([.5, .8, 1])
        clist = sc.gridcolors(3)
        clist = [clist[0], clist[1], clist[2]]
        colors = sc.objdict(treat50=clist[0], treat80=clist[1], treat100=clist[2])

        set_font(size=20)
        fig, axes = pl.subplots(2, 3, figsize=(20, 8))
        axes = axes.ravel()
        ai = 0

        # Plot symptomatic proportion
        ax = axes[ai]
        thisdf = df.loc[df['par'] == 'p_symp']
        thisdf['value'] = thisdf['value']*100
        sns.boxplot(data=thisdf, x="disease", y="value", hue="scenario", palette=clist, ax=ax)
        ax.set_title('Symptomatic share (F):\n% cases in women accompanied by VDS')
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_ylim(0, 100)
        # Turn off frame around legend
        ax.legend(frameon=False, prop={'size': 16})
        ai += 1

        # Plot symptomatic proportion
        ax = axes[ai]
        cs_df = sc.loadobj('results/zim_sti_care_seeking.obj')
        cs_df['p_symp_care'] = cs_df['p_symp_care']*100
        ax = sns.boxplot(data=cs_df, hue="scenario", y="p_symp_care", palette=clist, ax=ax)
        ax.set_title('Care-seeking share (F):\n% women with VDS who seek care')
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_ylim(0, 100)
        ax.get_legend().set_visible(False)
        ai += 1

        # Plot overall proportion treated
        ax = axes[ai]
        thisdf = df.loc[df['par'] == 'p_treat']
        thisdf['value'] = thisdf['value']*100
        ax = sns.boxplot(data=thisdf, x="disease", y="value", hue="scenario", palette=clist, ax=ax)
        ax.set_title('Treated share (F):\n% cases in women that are treated')
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_ylim(bottom=0)
        ax.get_legend().set_visible(False)
        ai += 1

        # Plot prevalence
        diseases = ['ng', 'ct', 'tv']
        t = results['treat50'][0]['time']
        si = sc.findfirst(t, 2010)
        ei = sc.findfirst(t, 2025)

        for i, disease in enumerate(diseases):
            ax = axes[i+3]

            for scenario in ['treat50', 'treat80', 'treat100']:
                res = results[scenario]
                prevalence = np.array([r[disease+'_prevalence'] for r in res])
                ax.plot(t[si:ei], prevalence.transpose()[si:ei, :]*100, alpha=0.6, lw=0.5, label=scenario, color=colors[scenario])

            ax.set_title(disease.upper()+' prevalence (%)')
            ax.set_ylim(bottom=0)
            ax.set_ylabel('')
            ax.set_xlabel('')

        fig.tight_layout()
        pl.savefig(f"figures/fig3_pars.png", dpi=100)

    print('Done.')
