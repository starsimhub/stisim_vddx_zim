"""
Create HIV model and interventions
"""

# %% Imports and settings
import numpy as np
import starsim as ss
import pandas as pd
import stisim as sti


def get_testing_products(location='zimbabwe'):
    """
    Define HIV products and testing interventions
    """
    # Load HIV test data:
    hiv_testing_data = pd.read_excel(f'data/{location}_20230725.xlsx', sheet_name='Testing & treatment', skiprows=1)
    HIV_tests_data_raw = hiv_testing_data.iloc[0:15, 1:43]
    HIV_low_cd4count_data_raw = hiv_testing_data.iloc[21:22, 1:43]
    HIV_low_cd4count_data_raw = HIV_low_cd4count_data_raw.iloc[:, 1:]
    HIV_tests_data_raw.index = HIV_tests_data_raw.iloc[:, 0]
    HIV_tests_data_raw = HIV_tests_data_raw.iloc[:, 1:]
    HIV_tests_data_raw.loc["Other_avg"] = HIV_tests_data_raw[HIV_tests_data_raw.index != "FSW"].mean()
    tivec = np.arange(start=1990, stop=2020 + 1, step=1)
    FSW_prop = np.interp(tivec,
                         HIV_tests_data_raw.loc["FSW"].index[~pd.isna(HIV_tests_data_raw.loc["FSW"].values)].astype(int),
                         HIV_tests_data_raw.loc["FSW"].values[~pd.isna(HIV_tests_data_raw.loc["FSW"].values)])
    other_prop = np.interp(tivec,
                           HIV_tests_data_raw.loc["Other_avg"].index[~pd.isna(HIV_tests_data_raw.loc["Other_avg"].values)].astype(int),
                           HIV_tests_data_raw.loc["Other_avg"].values[~pd.isna(HIV_tests_data_raw.loc["Other_avg"].values)])
    low_cd4count_prop = np.interp(tivec,
                                  HIV_low_cd4count_data_raw.iloc[0].index[~pd.isna(HIV_low_cd4count_data_raw.iloc[0].values)].astype(int),
                                  HIV_low_cd4count_data_raw.iloc[0].values[~pd.isna(HIV_low_cd4count_data_raw.iloc[0].values)])

    ####################################################################################################################
    # FSW Testing
    ####################################################################################################################

    # Eligible for testing are FSW agents who haven't been diagnosed or treated yet
    def fsw_eligibility(sim):
        return sim.networks.structuredsexual.fsw & ~sim.diseases.hiv.diagnosed & ~sim.diseases.hiv.on_art

    fsw_testing = sti.HIVTest(
        years=tivec,
        test_prob_data=FSW_prop,
        name='fsw_testing',
        eligibility=fsw_eligibility,
        label='fsw_testing',
    )

    ####################################################################################################################
    # Remaining population testing
    ####################################################################################################################

    # Eligible for testing are non-FSW agents who haven't been diagnosed or treated yet
    def other_eligibility(sim):
        return ~sim.networks.structuredsexual.fsw & ~sim.diseases.hiv.diagnosed & ~sim.diseases.hiv.on_art

    other_testing = sti.HIVTest(
        years=tivec,
        test_prob_data=other_prop,
        name='other_testing',
        eligibility=other_eligibility,
        label='other_testing',
    )

    ####################################################################################################################
    # Low CD4 count testing
    ####################################################################################################################
    # Eligible for testing are agents, who haven't been diagnosed yet and whose CD4 count is below 200.
    def low_cd4_eligibility(sim):
        return (sim.diseases.hiv.cd4 < 200) & ~sim.diseases.hiv.diagnosed

    low_cd4_testing = sti.HIVTest(
        years=tivec,
        test_prob_data=low_cd4count_prop,
        name='low_cd4_testing',
        eligibility=low_cd4_eligibility,
        label='low_cd4_testing',
    )

    return fsw_testing, other_testing, low_cd4_testing


def make_hiv():
    """ Make HIV arguments for sim"""
    hiv = sti.HIV(
        beta={'structuredsexual': [1, 1], 'maternal': [1, 0.]},
        beta_m2f=0.073,  # 0.101341,
        beta_f2m=0.025,  # 0.011625,
        beta_m2c=0.068,
        dur_on_art=ss.lognorm_ex(26, 5),
        init_prev_data=pd.read_csv('data/init_prev_hiv.csv'),
        rel_init_prev=0.5,
    )
    return [hiv]


def make_hiv_intvs(location='zimbabwe'):

    n_art = pd.read_csv(f'data/{location}_art.csv').set_index('year')
    n_vmmc = pd.read_csv(f'data/{location}_vmmc.csv').set_index('year')
    fsw_testing, other_testing, low_cd4_testing = get_testing_products()
    art = sti.ART(coverage_data=n_art)
    vmmc = sti.VMMC(coverage_data=n_vmmc)

    interventions=[
        fsw_testing,
        other_testing,
        low_cd4_testing,
        art,
        vmmc
    ]

    return interventions

