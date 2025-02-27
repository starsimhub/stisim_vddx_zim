"""
Plot degree of overtreatment by scenario
"""
import pandas as pd

# %% Imports and settings
import sciris as sc
import matplotlib.pyplot as pl
import seaborn as sns
from utils import set_font


if __name__ == '__main__':

    to_run = [
        # 'make_df',
        'make_plot',
    ]

    scenarios = ['treat50', 'treat80', 'treat100', 'treat50poc', 'treat80poc', 'treat100poc']
    scen_labels = {'treat50': 'Poor', 'treat80': 'Imperfect', 'treat100': 'Perfect', 'treat50poc': 'Poor (POC)', 'treat80poc': 'Imperfect (POC)', 'treat100poc': 'Perfect (POC)'}
    treatments = ['ng_tx', 'ct_tx', 'metronidazole']
    tx_labels = {'ng_tx':'NG', 'ct_tx':'CT', 'metronidazole':'MTNZ'}
    results = [tx+'.new_treated_unnecessary_f' for tx in treatments]
    results += [tx+'.new_treated_f' for tx in treatments]

    if 'make_df' in to_run:
        dfs = sc.autolist()
        for scenario in scenarios:
            sdf = sc.loadobj(f'results/{scenario}_sim.df')
            df = sdf.loc[:, sdf.columns.isin(['timevec']+results)]
            for tx in treatments:
                df[tx+'.overtx'] = df[tx+'.new_treated_unnecessary_f']/df[tx+'.new_treated_f']
            df['scenario'] = scen_labels[scenario]
            dfs += df
        bigdf = pd.concat(dfs)

        # Melt dataframe to long form
        dfm = bigdf.melt(id_vars=['timevec', 'scenario'], var_name='variable', value_name='value')
        dfm['treatment'] = dfm['variable'].apply(lambda x: x.split('.')[0])
        sc.saveobj(f'results/overtx.obj', dfm)

    if 'make_plot' in to_run:
        set_font(size=20)
        fig, axes = pl.subplots(2, 3, figsize=(20, 8))
        axes = axes.ravel()
        # clist = sc.vectocolor([0, .5, .8, 1], cmap='plasma_r')
        clist = sc.gridcolors(3)
        clist = [clist[0], clist[1], clist[2]]  #, clist[3]][1:]
        colors = sc.objdict(treat50=clist[0], treat80=clist[1], treat100=clist[2])

        # Care seeker management
        for pn, scenario in enumerate(['treat50', 'treat80', 'treat100']):
            ax = axes[pn]
            df = sc.loadobj(f'results/{scenario}_sim.df')
            dfplot = df.loc[(df.timevec >= 2010) & (df.timevec <= 2040)]
            x = dfplot.timevec
            sex = '_f'
            Y = [
                dfplot['syndromicmgmt.new_tx0'+sex],
                dfplot['syndromicmgmt.new_tx1'+sex],
                dfplot['syndromicmgmt.new_tx2'+sex],
                dfplot['syndromicmgmt.new_tx3'+sex],
            ]
            labels = ["0", "1", "2", "3"]
            ax.stackplot(x, *Y, baseline='zero', labels=labels, colors=sc.vectocolor(4, reverse=True))
            ax.set_title('Number of women presribed\nmultiple antibiotics')
            if pn==0: ax.legend(frameon=False, prop={'size': 12}, loc='upper left')
            ax.set_ylim(bottom=0)
            sc.SIticks(ax=ax)

        # Effect of POC
        df = sc.loadobj(f'results/overtx.obj')
        t = df.timevec.unique()
        si = sc.findfirst(t, 2010)
        ei = sc.findfirst(t, 2040)

        for pn, (txname, txlabel) in enumerate(tx_labels.items()):
            ax = axes[pn+3]
            for scenario in ['treat50', 'treat80', 'treat100']:
                socdf = df.loc[(df.scenario == scen_labels[scenario]) & (df.treatment == txname) & (df.variable == txname+'.new_treated_unnecessary_f')]
                socy = socdf['value'][si:ei]
                socy = socy.rolling(10, min_periods=1).mean()
                ax.plot(t[si:ei], socy, label=scen_labels[scenario], color=colors[scenario])
            for scenario in ['treat50', 'treat80', 'treat100']:
                pocdf = df.loc[(df.scenario == scen_labels[scenario+'poc']) & (df.treatment == txname) & (df.variable == txname+'.new_treated_unnecessary_f')]
                pocy = pocdf['value'][si:ei]
                pocy = pocy.rolling(10, min_periods=1).mean()
                ax.plot(t[si:ei], pocy, label=scen_labels[scenario], color=colors[scenario], ls='--')

            if pn == 0:
                h,l = ax.get_legend_handles_labels()  # #Get the legend handles and lables
                l1 = ax.legend(h[:3], l[:3], loc='upper left', frameon=False, prop={'size': 12})
                from matplotlib.lines import Line2D
                myHandle = [Line2D([], [], ls='-', color='k'), Line2D([], [], ls='--', color='k')]
                l2 = ax.legend(handles=myHandle, labels=['SOC', 'POC'], loc='upper left', bbox_to_anchor=(0.3, 1), frameon=False, prop={'size': 12})
                ax.add_artist(l1)

            ax.set_title(f'{txlabel} overtreatment')
            ax.set_ylim(bottom=0)
            sc.SIticks(ax)

        fig.tight_layout()
        pl.savefig(f"figures/fig4_overtx.png", dpi=100)
 
    print('Done!')





