# Script for generating main text & supplementary figures

import numpy as np
import pandas as pd
import seaborn as sb
# import seaborn.objects as so
from matplotlib import pyplot as plt
from matplotlib.transforms import ScaledTranslation
import os, sys

# ============================================================================ #
# return high 6 bits
def agent_state(id):
    return id >> np.uint64(58)
# ============================================================================ #
# zero out the high 6 bits
def agent_id(id):
    return id & ~(np.uint64(0b00111111) << np.uint64(58))
# ============================================================================ #
def column_names(file, length):
    tmp = file.read(length)
    return tmp[0:-1].decode("ascii").split("\x00")
# ============================================================================ #
def select_transitions(data, disease_state, context_lo=-1, context_hi=255):
    return data[(data["state"] == disease_state) &
        ((context_lo < data["context"]) & (data["context"] <= context_hi))]
# ============================================================================ #
def read_event_file(filename):
    event_t = np.dtype({
        "names": ["agent_id","location_id","timestep","context","state",
            "variant"],
        "formats": [np.uint64, np.uint64, np.uint16, np.uint8, np.uint8,
            np.uint8]
    }, align=True)

    with open(filename, mode="rb") as file:
        hdr = np.frombuffer(file.read(48), dtype=np.uint64, count=6)
        demog_names = column_names(file, hdr[4])
        col_names = column_names(file, hdr[5])

        col_names = ["agent_id", "location_id", "timestep",
            "infection_context", "disease_state", "variant"]

        tract_fips = np.frombuffer(file.read(np.uint64(8) * hdr[0]),
            dtype=np.uint64, count=hdr[0])
    
        demog_table = np.frombuffer(file.read(np.uint64(4) * hdr[0] * hdr[3]),
            dtype=np.uint32, count=hdr[0] * hdr[3])

        pos = file.tell()
        file.seek(0, 2)
        n_event = np.uint64((file.tell() - pos) / event_t.itemsize)
        file.seek(pos, 0)

        data = np.frombuffer(file.read(n_event * np.uint64(event_t.itemsize)),
            dtype=event_t, count=n_event)

        return demog_names, demog_table, col_names, data
# ============================================================================ #
# Function to generate figure 4 - household size effects
# events_files_all: list of events file IDs ('run_0{ID}.events.bin')
# state: state FIPS code (e.g. New Mexico - 35)
# def figure_4_gen(events_files_all,state):
def figure_4_gen(data_dir, state=35, compare="race", SA=False):

    # UrbanPop population analysis
    if state == 35:
        agent_file = os.path.join(data_dir, "aux", "nm.csv")
        events_files_all = list(range(1,11))
    elif state == 12:
        agent_file = os.path.join(data_dir, "aux", "fl.csv") # supplement figure for FL
        events_files_all = list(range(21,31))
    else:
        raise Exception('State file not available - only 35 (NM) and 12 (FL) are supported')

    if not compare in ["race","ethnicity"]:
        raise Exception("Invalid compare variable: MUST be either \"race\" or \"ethnicity\"")

    data_all =  pd.read_csv(agent_file, index_col= 'person_id', 
                            dtype={'fips_code' : np.uint64, 
                                   'person_id' : np.uint32, 
                                   'household_id' : np.uint32, 
                                   'person_race' : np.uint8,
                                   'person_ethnicity' : np.uint8, 
                                   'household_size' : np.uint8})

    aian_total = np.count_nonzero(data_all['person_race'] == 3)
    white_total = np.count_nonzero(data_all['person_race'] == 0)
    hisp_total = np.count_nonzero(data_all['person_ethnicity'] == 1)
    non_hisp_total = np.count_nonzero(data_all['person_ethnicity'] == 0)

    HH_total = [0] * 7
    for h in range(1,7):
        HH_total[h-1] = np.count_nonzero(data_all['household_size'] == h)
    HH_total[6] = np.count_nonzero(data_all['household_size'] > 6)

    data_vis = pd.DataFrame()

    total_fips = data_all['fips_code']
    total_fips = np.unique(total_fips,return_counts=True)

    data_vis['FIPS'] = total_fips[0]
    data_vis['Total'] = total_fips[1]

    hh_fips = []

    for index, row in data_all.iterrows():
        hh_fips.extend([row['fips_code']] * row['household_size'])
    #     for g in range(row['household_size']):
    #         hh_fips.append(row['fips_code'])

    hh_fips = np.unique(hh_fips,return_counts=True)
    hh_fips_dict = dict(zip(hh_fips[0],hh_fips[1]))
    hh_fips_tmp = np.empty(len(total_fips[0]),int)

    for f in range(0,len(total_fips[0])):
        hh_fips_tmp[f] = hh_fips_dict[total_fips[0][f]]

    data_vis['HH size'] = hh_fips_tmp


    # EpiCast simulation analysis
    exposed_state = 1

    exposed_white = []
    exposed_aian = []

    exposed_hisp = []
    exposed_non_hisp = []

    exposed_white_overall = []
    exposed_aian_overall = []

    exposed_hisp_overall = []
    exposed_non_hisp_overall = []

    exposed_HH = []

    for e in events_files_all:
        if not SA:
            events_file = os.path.join(data_dir, "main_figures", f"run_{e:03}.events.bin")
        else:
            events_file = os.path.join(data_dir, "SA", f"run_{e:03}.events.bin")

        dn, dt, cn, events = read_event_file(events_file)

        # Select household exposure only
        exposed = select_transitions(events, exposed_state, context_lo=-1, context_hi=0) 

        full_ids = exposed['agent_id']
        all_states = agent_state(full_ids)

        ids = agent_id(full_ids[all_states == state])

        data = data_all.loc[ids]

        exposed_non_hisp_temp = np.count_nonzero(data['person_ethnicity'] == 0)
        exposed_hisp_temp = np.count_nonzero(data['person_ethnicity'] == 1)

        exposed_aian_temp = np.count_nonzero(data['person_race'] == 3)
        exposed_white_temp = np.count_nonzero(data['person_race'] == 0)

        exposed_hisp.append(exposed_hisp_temp / hisp_total)
        exposed_non_hisp.append(exposed_non_hisp_temp / non_hisp_total)

        exposed_aian.append(exposed_aian_temp / aian_total)
        exposed_white.append(exposed_white_temp / white_total)

        # Select overall exposure
        exposed = select_transitions(events, exposed_state, context_lo=-1, context_hi=255) 

        full_ids = exposed['agent_id']
        all_states = agent_state(full_ids)

        ids = agent_id(full_ids[all_states == state])

        data = data_all.loc[ids]

        exposed_non_hisp_temp = np.count_nonzero(data['person_ethnicity'] == 0)
        exposed_hisp_temp = np.count_nonzero(data['person_ethnicity'] == 1)

        exposed_aian_temp = np.count_nonzero(data['person_race'] == 3)
        exposed_white_temp = np.count_nonzero(data['person_race'] == 0)

        for h in range(1,7):
            exposed_HH.append(np.count_nonzero(data['household_size'] == h) / HH_total[h-1])
        exposed_HH.append(np.count_nonzero(data['household_size'] > 6) / HH_total[6])

        exposed_hisp_overall.append(exposed_hisp_temp / hisp_total)
        exposed_non_hisp_overall.append(exposed_non_hisp_temp / non_hisp_total)

        exposed_aian_overall.append(exposed_aian_temp / aian_total)
        exposed_white_overall.append(exposed_white_temp / white_total)

        exposed_fips = data['fips_code']
        exposed_fips = np.unique(exposed_fips,return_counts=True)

        exposed_dict = dict(zip(exposed_fips[0],exposed_fips[1]))

        exposed_fips_tmp = np.empty(len(total_fips[0]), dtype=int)

        for f in range(0,len(total_fips[0])):
            if total_fips[0][f] in exposed_fips[0]:
                exposed_fips_tmp[f] = exposed_dict[total_fips[0][f]]
            else:
                exposed_fips_tmp[f] = 0
    
        data_vis['Exposed' + str(e)] = exposed_fips_tmp

    # Aggregate to census tract level
    data_vis = data_vis.groupby(by=["FIPS"]).sum()

    filter_col = [col for col in data_vis if col.startswith('Exposed')]
    data_vis['Exposed'] = data_vis[filter_col].apply(np.mean,axis=1)

    data_vis['%Exposed'] = data_vis['Exposed'] / data_vis['Total']
    data_vis['Mean Household Size'] = data_vis['HH size'] / data_vis['Total']


    plt.rcParams.update({'font.size': 16})

    fig, axes = plt.subplots(nrows=1, ncols=3)

    fig.set_figheight(3)
    fig.set_figwidth(12)

    #Figure 4A
    hh_size_dens_race = pd.DataFrame()
    hh_size_dens_race['Race'] = data_all['person_race']
    race_map = {0:'White', 1:'Black', 2:'Asian', 3:'AIAN', 4:'NHPI',5:'Other',6:'Multi'}
    hh_size_dens_race['Race'] = [race_map[h] for h in hh_size_dens_race['Race']]
    hh_size_dens_race['Household (HH) Size'] = data_all['household_size']
    hh_size_dens_race = hh_size_dens_race[(hh_size_dens_race['Race'] == 'White') | 
                                          (hh_size_dens_race['Race'] == 'AIAN')]

    hh_size_dens_eth = pd.DataFrame()
    hh_size_dens_eth['Ethnicity'] = data_all['person_ethnicity']
    ethnicity_map = {0:'Non-Hisp.', 1:'Hisp.'}
    hh_size_dens_eth['Ethnicity'] = [ethnicity_map[h] for h in hh_size_dens_eth['Ethnicity']]
    hh_size_dens_eth['Household (HH) Size'] = data_all['household_size']

    cmap = sb.color_palette("blend:#F8766D,#00BFC4", n_colors=2)

    if compare == "race":
        hue_order = ['AIAN','White']
        sb.kdeplot(data=hh_size_dens_race, x='Household (HH) Size', hue='Race', 
                hue_order=hue_order, common_norm=False, bw_adjust=5.5, clip=(1,15), 
                fill=True, alpha=0.2,  palette=cmap, ax=axes[0])
    elif compare == "ethnicity":
        hue_order = ['Non-Hisp.','Hisp.']
        sb.kdeplot(data=hh_size_dens_eth, x='Household (HH) Size', hue='Ethnicity', 
            common_norm=False, bw_adjust=8, clip=(1,15), 
            fill=True, alpha=0.2, palette=cmap, ax=axes[0])

    axes[0].set_xlim(left=1)

    axes[0].text(
        0.0, 1.0, 'A)', transform=(
            axes[0].transAxes + ScaledTranslation(-20/72, +7/72, fig.dpi_scale_trans)),
        fontsize='medium', va='bottom', fontfamily='sans-serif')

    #Figure 4B

    # sb.scatterplot(data=data_vis, x="Mean Household Size", y="%Exposed", alpha=0.5, ax=axes[1])
    # axes[1].set_xlabel('Mean household size among\n individuals living in tract')
    # # axes[1].set_ylabel('Proportion infected\n w/in household')
    # axes[1].set_ylabel('Proportion infected\n overall')

    axes[1].text(
        0.0, 1.0, 'B)', transform=(
            axes[1].transAxes + ScaledTranslation(-20/72, +7/72, fig.dpi_scale_trans)),
        fontsize='medium', va='bottom', fontfamily='sans-serif')

    hh_size = pd.DataFrame()
    hh_size['Household (HH) Size'] = ([str(h) for h in list(range(1,7))] + ['7+']) * len(exposed_white)
    hh_size['Exposure'] = exposed_HH

    sb.barplot(hh_size,x='Household (HH) Size',y='Exposure',errorbar=None, ax=axes[1])

    sb.stripplot(hh_size,x='Household (HH) Size',y='Exposure', alpha=0.5, s=5.5, ax=axes[1], 
                 color='black')
    
    axes[1].set_ylabel('Proportion infected')
    
    #Figure 4C
    palette = ['#fde3e1','#F8766D','#ccf2f3','#00BFC4']

    HH_exp_groups = pd.DataFrame()
    if compare == "race":
        HH_exp_groups['Racial/ethnic group'] = ['White']*len(exposed_white)*2 + ['AIAN']*len(exposed_white)*2
        HH_exp_groups['Context'] = ['Non-Household']*len(exposed_white) + ['Household']*len(exposed_white) + ['Non-Household']*len(exposed_white) + ['Household']*len(exposed_white)
        HH_exp_groups['Exposure'] =  exposed_white_overall + exposed_white + exposed_aian_overall + exposed_aian
    elif compare == "ethnicity":
        HH_exp_groups['Racial/ethnic group'] = ['Non-Hisp.']*len(exposed_white)*2 + ['Hisp.']*len(exposed_white)*2
        HH_exp_groups['Context'] = ['Non-Household']*len(exposed_white) + ['Household']*len(exposed_white) + ['Non-Household']*len(exposed_white) + ['Household']*len(exposed_white)
        HH_exp_groups['Exposure'] = exposed_non_hisp_overall + exposed_non_hisp + exposed_hisp_overall + exposed_hisp

    sb.barplot(HH_exp_groups,x='Racial/ethnic group',y='Exposure',errorbar=None, ax=axes[2], 
               palette=palette, hue='Context', order=hue_order, dodge=False)
    sb.stripplot(HH_exp_groups,x='Racial/ethnic group',y='Exposure', alpha=0.5, s=7, ax=axes[2], 
                 color='black', order=hue_order)
    
    for bars, colors in zip(axes[2].containers, (palette[0::2], palette[1::2])):
     for bar, color in zip(bars, colors):
          bar.set_facecolor(color)

    axes[2].legend_.remove()

    # print(HH_exp_groups.groupby(by=['Racial/ethnic group','Context']).mean())
    # print(HH_exp_groups.groupby(by=['Racial/ethnic group','Context']).var())

    # axes[2].text(
    #     0.0, 1.0, 'C)', transform=(
    #         axes[2].transAxes + ScaledTranslation(-20/72, +7/72, fig.dpi_scale_trans)),
    #     fontsize='medium', va='bottom', fontfamily='sans-serif')
    
    # axes[2].set_ylabel('Proportion infected\n w/in household')
    axes[2].set_ylabel('Proportion infected')

    labels = [['Non-HH','Non-HH'], ['HH','HH']]

    # for c in range(len(axes[2].containers)):
    axes[2].bar_label(axes[2].containers[0], labels=labels[0], padding=-30)
    axes[2].bar_label(axes[2].containers[1], labels=labels[1], label_type='center') # label_type='center',

    axes[2].text(
        0.0, 1.0, 'C)', transform=(
            axes[2].transAxes + ScaledTranslation(-20/72, +7/72, fig.dpi_scale_trans)),
        fontsize='medium', va='bottom', fontfamily='sans-serif')

    axes[0].spines[['right','top']].set_visible(False)
    axes[1].spines[['right','top']].set_visible(False)
    axes[2].spines[['right','top']].set_visible(False)

    plt.tight_layout()
    if not SA:
        ofile = "fig4_" + str(state) + "_" + compare + ".pdf"
    else:
        ofile = "fig4_" + str(state) + "_" + compare + "_SA.pdf"
    plt.savefig(ofile)

    return
# ============================================================================ #
# Function to generate figure 5 - workplace effects
# events_files_all: dictionary of events file IDs ('run_0{ID}.events.bin') 
#   as keys to model numbers as values (e.g. {7:3,8:3,9:2,10:2})
# state: state FIPS code (e.g. New Mexico - 35)
# def figure_5_gen(events_files_all,state):
def figure_5_gen(data_dir, SA=False):

    state = 36

    high_risk_industries = [445, 446, 448, 451, 452, 622, 623, 624, 713, 722]

    # UrbanPop population analysis
    agent_file = os.path.join(data_dir, "aux", "ny.csv")
    data_all =  pd.read_csv(agent_file, index_col= 'person_id', 
                            dtype={'fips_code' : np.uint64, 
                                   'person_id' : np.uint32, 
                                   'person_race' : np.uint8, 
                                   'person_naics' : np.uint16})

    events_files_all = {k : 1 for k in range(11, 21)}
    events_files_all.update({k : 3 for k in range(31, 41)})

    race_total = data_all['person_race']
    race_total = np.unique(race_total,return_counts=True)
    race_total = dict(zip(race_total[0],race_total[1]))

    high_risk = data_all['person_naics']
    high_risk = np.unique(high_risk,return_counts=True)
    high_risk = dict(zip(high_risk[0],high_risk[1]))

    high_risk_total = 0
    low_risk_total = 0

    for i in high_risk:
        if i in high_risk_industries:
            high_risk_total += high_risk[i]
        else:
            low_risk_total += high_risk[i]


    black_total = race_total[1]
    white_total = race_total[0] 

    race_naics_total = zip(data_all['person_race'], data_all['person_naics'])
    race_naics_total = np.unique(list(race_naics_total),return_counts=True,axis=0)
    race_naics_total = dict(zip([tuple(h) for h in race_naics_total[0]],race_naics_total[1]))

    high_risk_black = 0
    high_risk_white = 0
    
    employed_black = 0
    employed_white = 0
  
    for i in race_naics_total:
        if i[0] == 0 and i[1] in high_risk_industries:
            high_risk_white += race_naics_total[i]
        if i[0] == 1 and i[1] in high_risk_industries:
            high_risk_black += race_naics_total[i]
        if i[0] == 0 and i[1] > 1:
            employed_white += race_naics_total[i]
        if i[0] == 1 and i[1] > 1:
            employed_black += race_naics_total[i]

    # Overall black/white population as denominator
    # high_risk_black = high_risk_black / black_total
    # high_risk_white = high_risk_white / white_total

    # Employed black/white population as denominator
    high_risk_black = high_risk_black / employed_black
    high_risk_white = high_risk_white / employed_white

    data_vis = pd.DataFrame()

    total_fips = data_all['fips_code']
    total_fips = np.unique(total_fips,return_counts=True)

    data_vis['FIPS'] = total_fips[0]
    data_vis['Total'] = total_fips[1]

    fips_naics_total = zip(data_all['fips_code'],data_all['person_naics'])
    fips_naics_total = np.unique(list(fips_naics_total),return_counts=True,axis=0)
    fips_naics_total = dict(zip([tuple(h) for h in fips_naics_total[0]],fips_naics_total[1]))

    high_risk_fips_tmp = np.empty(len(total_fips[0]),int)
    industry_fips_tmp = np.empty(len(total_fips[0]),int)

    high_risk_fips_dict = {}
    industry_fips_dict = {}

    for i in fips_naics_total:
        if i[1] in high_risk_industries:
            if i[0] in high_risk_fips_dict:
                high_risk_fips_dict[i[0]] += fips_naics_total[i]
            else:
                high_risk_fips_dict[i[0]] = fips_naics_total[i]
        if i[1] > 1:
            if i[0] in industry_fips_dict:
                industry_fips_dict[i[0]] += fips_naics_total[i]
            else:
                industry_fips_dict[i[0]] = fips_naics_total[i]

    for f in range(0,len(total_fips[0])):
        if total_fips[0][f] in high_risk_fips_dict:
            high_risk_fips_tmp[f] = high_risk_fips_dict[total_fips[0][f]]
        else:
            high_risk_fips_tmp[f] = 0
        if total_fips[0][f] in industry_fips_dict:
            industry_fips_tmp[f] = industry_fips_dict[total_fips[0][f]]
        else:
            industry_fips_tmp[f] = 0

    data_vis['Individuals living in tract employed in a high risk industry'] = high_risk_fips_tmp
    data_vis['Individuals living in tract employed'] = industry_fips_tmp

    # EpiCast simulation analysis
    exposed_state = 1

    exposed_white = []
    exposed_black = []

    exposed_high_risk = []
    exposed_low_risk = []

    exposed_white_overall_1 = []
    exposed_black_overall_1 = []
    exposed_white_overall_2 = []
    exposed_black_overall_2 = []


    for e in events_files_all:
        if not SA:
            events_file = os.path.join(data_dir, "main_figures", f"run_{e:03}.events.bin")
        else:
            events_file = os.path.join(data_dir, "SA", f"run_{e:03}.events.bin")
        

        dn, dt, cn, events = read_event_file(events_file)
        if events_files_all[e] == 1: # check if model 1

            # Select workplace exposure only
            exposed = select_transitions(events, exposed_state, context_lo=3, context_hi=4) 

            full_ids = exposed['agent_id']
            all_states = agent_state(full_ids)

            # data = np.empty(len(full_ids), dtype=get_agent_type())

            ids = agent_id(full_ids[all_states == state])

            data = data_all.loc[ids]

            exposed_black_temp = np.count_nonzero(data['person_race'] == 1)
            exposed_white_temp = np.count_nonzero(data['person_race'] == 0)

            # Overall black/white population as denominator
            # exposed_black.append(exposed_black_temp / black_total)
            # exposed_white.append(exposed_white_temp / white_total)

            # Employed black/white population as denominator
            exposed_black.append(exposed_black_temp / employed_black)
            exposed_white.append(exposed_white_temp / employed_white)

            # exposed_high_risk_temp = 0
            # exposed_low_risk_temp = 0

            # for i in high_risk:
            #     if i in high_risk_industries:
            #         exposed_high_risk_temp += np.count_nonzero(data['person_naics'] == i)
            #     else:
            #         exposed_low_risk_temp += np.count_nonzero(data['person_naics'] == i)

            # exposed_high_risk.append(exposed_high_risk_temp / high_risk_total)
            # exposed_low_risk.append(exposed_low_risk_temp / low_risk_total)

            # Measure workplace exposure in each census tract
            # exposed_fips = data['fips_code']
            # exposed_fips = np.unique(exposed_fips,return_counts=True)

            # exposed_dict = dict(zip(exposed_fips[0],exposed_fips[1]))

            # exposed_fips_tmp = np.empty(len(total_fips[0]), dtype=int)

            # for f in range(0,len(total_fips[0])):
            #     if total_fips[0][f] in exposed_fips[0]:
            #         exposed_fips_tmp[f] = exposed_dict[total_fips[0][f]]
            #     else:
            #         exposed_fips_tmp[f] = 0
        
            # data_vis['Exposed' + str(e)] = exposed_fips_tmp

            # Select all exposure events (regardless of context)       
            exposed_overall = select_transitions(events, exposed_state, context_lo=-1, context_hi=255) 

            full_ids = exposed_overall['agent_id']
            all_states = agent_state(full_ids)

            ids = agent_id(full_ids[all_states == state])

            data = data_all.loc[ids]

            exposed_black_temp = np.count_nonzero(data['person_race'] == 1)
            exposed_white_temp = np.count_nonzero(data['person_race'] == 0)

            exposed_black_overall_1.append(exposed_black_temp / black_total)
            exposed_white_overall_1.append(exposed_white_temp / white_total)

            # # Measure overall exposure in each census tract
            exposed_fips = data['fips_code']
            exposed_fips = np.unique(exposed_fips,return_counts=True)

            exposed_dict = dict(zip(exposed_fips[0],exposed_fips[1]))

            exposed_fips_tmp = np.empty(len(total_fips[0]), dtype=int)

            for f in range(0,len(total_fips[0])):
                if total_fips[0][f] in exposed_fips[0]:
                    exposed_fips_tmp[f] = exposed_dict[total_fips[0][f]]
                else:
                    exposed_fips_tmp[f] = 0
        
            data_vis['Exposed' + str(e)] = exposed_fips_tmp

            exposed_high_risk_temp = 0
            exposed_low_risk_temp = 0

            for i in high_risk:
                if i in high_risk_industries:
                    exposed_high_risk_temp += np.count_nonzero(data['person_naics'] == i)
                else:
                    exposed_low_risk_temp += np.count_nonzero(data['person_naics'] == i)

            exposed_high_risk.append(exposed_high_risk_temp / high_risk_total)
            exposed_low_risk.append(exposed_low_risk_temp / low_risk_total)

        else: # model 2

            # Select all exposure events (regardless of context)
            exposed_overall = select_transitions(events, exposed_state, context_lo=-1, context_hi=255)

            full_ids = exposed_overall['agent_id']
            all_states = agent_state(full_ids)

            ids = agent_id(full_ids[all_states == state])

            data = data_all.loc[ids]

            exposed_black_temp = np.count_nonzero(data['person_race'] == 1)
            exposed_white_temp = np.count_nonzero(data['person_race'] == 0)

            exposed_black_overall_2.append(exposed_black_temp / black_total)
            exposed_white_overall_2.append(exposed_white_temp / white_total)

    # Aggregate to census tract level
    data_vis = data_vis.groupby(by=["FIPS"]).sum()

    filter_col = [col for col in data_vis if col.startswith('Exposed')]
    data_vis['Exposed'] = data_vis[filter_col].apply(np.mean,axis=1)

    # Denominator: all individuals in tract
    data_vis['%Exposed'] = data_vis['Exposed'] / data_vis['Total']

    # Denominator: all workers in tract
    # data_vis['%Exposed'] = data_vis['Exposed'] / data_vis['Individuals living in tract employed']

    # Denominator: all individuals in tract
    data_vis['Proportion of individuals living in tract working in a high risk industry'] =  data_vis['Individuals living in tract employed in a high risk industry'] / data_vis['Total']
    
    # Denominator: all workers in tract
    # data_vis['Proportion of individuals living in tract working in a high risk industry'] =  data_vis['Individuals living in tract employed in a high risk industry'] / data_vis['Individuals living in tract employed']
    
    plt.rcParams.update({'font.size': 16})

    fig, axes = plt.subplots(nrows=2, ncols=2)
    fig.set_figheight(11.25)
    fig.set_figwidth(15)

    cmap = sb.color_palette("blend:#F8766D,#00BFC4", n_colors=2)

    # Figure 5A
    high_risk_race = pd.DataFrame()
    high_risk_race['Race'] = ['Black', 'White']
    high_risk_race['Proportion employed in High Risk Industries'] = [high_risk_black, high_risk_white]

    sb.barplot(high_risk_race,x='Race',y='Proportion employed in High Risk Industries',errorbar=None, palette=cmap, ax=axes[0][0])
    axes[0][0].set_xlabel('Racial/ethnic group')
    
    # Denominator: all individuals in racial/ethnic group
    # axes[0][0].set_ylabel('Proportion employed in \nhigh risk industries')
    
    # Denominator: all workers in racial/ethnic group
    axes[0][0].set_ylabel('Proportion of workers employed in \nhigh risk industries')

    axes[0][0].text(
        0.0, 1.0, 'A)', transform=(
            axes[0][0].transAxes + ScaledTranslation(-20/72, +7/72, fig.dpi_scale_trans)),
        fontsize='medium', va='bottom', fontfamily='sans-serif')

    # Figure 5B

    sb.scatterplot(data=data_vis, 
                   x="Proportion of individuals living in tract working in a high risk industry", 
                   y="%Exposed", alpha=0.5, ax=axes[0][1])
    
    # # Portion of individuals in census tract infected in workplace
    # axes[0][1].set_xlim(right=0.35)
    # axes[0][1].set_xlabel('Proportion living in tract employed in high risk industries')
    # axes[0][1].set_ylabel('Proportion infected w/in workplace')

    # # Portion of individuals in census tract infected
    axes[0][1].set_xlim(right=0.35)
    axes[0][1].set_xlabel('Proportion living in tract employed in high risk industries')
    axes[0][1].set_ylabel('Proportion infected overall')

    # Portion of infected workers in census tract
    # axes[0][1].set_xlim(right=0.6)
    # axes[0][1].set_xlabel('Proportion of workers living in tract employed \nin high risk industries')
    # axes[0][1].set_ylabel('Proportion infected overall')

    # Portion of workers infected at work in census tract
    # axes[0][1].set_xlim(right=0.6)
    # axes[0][1].set_xlabel('Proportion of workers living in tract employed \nin high risk industries')
    # axes[0][1].set_ylabel('Proportion of workers infected \nw/in workplace')

    axes[0][1].text(
        0.0, 1.0, 'B)', transform=(
            axes[0][1].transAxes + ScaledTranslation(-20/72, +7/72, fig.dpi_scale_trans)),
        fontsize='medium', va='bottom', fontfamily='sans-serif')

    palette = ['#e69138','#cc0000']

    # high_risk_exp_groups = pd.DataFrame()
    # high_risk_exp_groups['Work type'] = ['High Risk']*len(exposed_white) + ['Other']*len(exposed_white)
    # high_risk_exp_groups['Exposure'] = exposed_high_risk + exposed_low_risk

    # sb.barplot(high_risk_exp_groups,x='Work type',
    #            y='Exposure',errorbar=None, ax=axes[0][1],palette=palette)
    # sb.stripplot(high_risk_exp_groups,x='Work type',
    #              y='Exposure', alpha=0.5, s=10, color='black', ax=axes[0][1])
    
    # axes[0][1].set_xlabel('Employment type')

    # Denominator: all individuals in racial/ethnic group
    # axes[1][0].set_ylabel('Proportion infected w/in workplace')

    # Denominator: all workers in racial/ethnic group
    # axes[0][1].set_ylabel('Proportion of workers infected \nw/in workplace')
    # axes[0][1].set_ylabel('Proportion of workers infected')

    # axes[0][1].text(
    #     0.0, 1.0, 'B)', transform=(
    #         axes[0][1].transAxes + ScaledTranslation(-20/72, +7/72, fig.dpi_scale_trans)),
    #     fontsize='medium', va='bottom', fontfamily='sans-serif')

    # Figure 5C

    high_risk_exp_groups = pd.DataFrame()
    high_risk_exp_groups['Racial/ethnic group'] = ['Black']*len(exposed_white) + ['White']*len(exposed_white)
    high_risk_exp_groups['Exposure'] = exposed_black + exposed_white

    sb.barplot(high_risk_exp_groups,x='Racial/ethnic group',
               y='Exposure',errorbar=None, palette=cmap, ax=axes[1][0])
    sb.stripplot(high_risk_exp_groups,x='Racial/ethnic group',
                 y='Exposure', alpha=0.5, s=10, color='black', ax=axes[1][0])
    
    axes[1][0].set_xlabel('Racial/ethnic group')

    # Denominator: all individuals in racial/ethnic group
    # axes[1][0].set_ylabel('Proportion infected w/in workplace')

    # Denominator: all workers in racial/ethnic group
    axes[1][0].set_ylabel('Proportion of workers infected \nw/in workplace')

    axes[1][0].text(
        0.0, 1.0, 'C)', transform=(
            axes[1][0].transAxes + ScaledTranslation(-20/72, +7/72, fig.dpi_scale_trans)),
        fontsize='medium', va='bottom', fontfamily='sans-serif')


    # Figure 5D

    # Overall infections under each model
    high_risk_exp_groups = pd.DataFrame()
    high_risk_exp_groups['Racial/ethnic group'] = ['Black']*len(exposed_white_overall_1) + ['White']*len(exposed_white_overall_1) + ['Black']*len(exposed_white_overall_2) + ['White']*len(exposed_white_overall_2)
    high_risk_exp_groups['Model'] = ['Model 1']*len(exposed_white_overall_1)*2 + ['Model 2']*len(exposed_white_overall_1)*2
    high_risk_exp_groups['Exposure'] = exposed_black_overall_1 + exposed_white_overall_1 + exposed_black_overall_2 + exposed_white_overall_2

    # print(high_risk_exp_groups.groupby(by=['Racial/ethnic group','Model']).mean())
    # print(high_risk_exp_groups.groupby(by=['Racial/ethnic group','Model']).var())

    # sb.barplot(high_risk_exp_groups,x='Racial/ethnic group',
    #            y='Exposure', hue='Model', errorbar=None, 
    #            palette="YlGnBu", ax=axes[1][1])
    # sb.stripplot(high_risk_exp_groups,x='Racial/ethnic group',
    #              y='Exposure', hue='Model', alpha=0.5, s=10, 
    #              color='black', ax=axes[1][1], dodge=True)

    # axes[1][1].set_xlabel('Racial/ethnic group')
    # axes[1][1].set_ylabel('Proportion infected overall')

    # Overall difference in infection rate between black & white pop.
    high_risk_exp_groups = pd.DataFrame()
    high_risk_exp_groups['Model'] = ['Model 1']*len(exposed_white_overall_1) + ['Model 2']*len(exposed_white_overall_1)
    high_risk_exp_groups['Exposure'] =  [ exposed_black_overall_1[h] - exposed_white_overall_1[h] for h in range(len(exposed_white_overall_1))] + [ exposed_black_overall_2[h] - exposed_white_overall_2[h] for h in range(len(exposed_white_overall_2))]

    sb.barplot(high_risk_exp_groups,x='Model',
               y='Exposure', errorbar=None, 
               palette="YlGnBu", ax=axes[1][1])
    sb.stripplot(high_risk_exp_groups,x='Model',
                 y='Exposure', alpha=0.5, s=10, 
                 color='black', ax=axes[1][1], dodge=True)

    axes[1][1].set_xlabel('Model')
    # axes[1][1].set_ylabel('Overall difference (%Black pop. infected - \n%White pop. infected)')
    axes[1][1].set_ylabel('Δ proportion infected overall\n(Black - White)')


    # handles, labels = axes[1][1].get_legend_handles_labels()
    # axes[1][1].legend(handles[0:2], labels[0:2])

    axes[1][1].text(
        0.0, 1.0, 'D)', transform=(
            axes[1][1].transAxes + ScaledTranslation(-20/72, +7/72, fig.dpi_scale_trans)),
        fontsize='medium', va='bottom', fontfamily='sans-serif')

    
    axes[0][0].spines[['right','top']].set_visible(False)
    axes[0][1].spines[['right','top']].set_visible(False)
    axes[1][0].spines[['right','top']].set_visible(False)
    axes[1][1].spines[['right','top']].set_visible(False)

    plt.tight_layout()
    if not SA:
        plt.savefig('fig5.pdf')
    else:
        plt.savefig('figS5.pdf')

    return

# ============================================================================ #

# Function to generate figure 6 - school effects
# events_files_all: dictionary of events file IDs ('run_0{ID}.events.bin') as keys to mode$
# state: state FIPS code (e.g. New Mexico - 35)
# def figure_6_gen(events_files_all,state):
def figure_6_gen(data_dir, supplement=False, SA=False):

    state = 12

    # UrbanPop population analysis
    agent_file = os.path.join(data_dir, "aux", "fl.csv")
    data_all =  pd.read_csv(agent_file, index_col= 'person_id', 
                            dtype={'fips_code' : np.uint64, 'person_id' : np.uint32, 
                                'household_id' : np.uint32, 'person_race' : np.uint8,
                                'person_ethnicity' : np.uint8, 
                                'household_size' : np.uint8})

    events_files_all = {k : 1 for k in range(21, 31)}
    events_files_all.update({k : 3 for k in range(41, 51)})

    eth_total = data_all['person_ethnicity']
    eth_total = np.unique(eth_total,return_counts=True)
    eth_total = dict(zip(eth_total[0],eth_total[1]))

    hispanic_total = eth_total[1]
    non_hispanic_total = eth_total[0]

    fl_schoolgroups = pd.read_csv(os.path.join(data_dir, "aux", "fl_agent_schoolgroups.csv"))

    fl_school_eth = np.unique(fl_schoolgroups['ethnicity'], return_counts=True)
    
    # Total hispanic/non-hispanic school attendees
    hispanic_school_total = fl_school_eth[1][1]
    non_hispanic_school_total = fl_school_eth[1][0]

    fl_school = np.unique(fl_schoolgroups['census_tracts'], return_counts=True)
    fl_school = dict(zip(fl_school[0],fl_school[1]))

    data_vis = pd.DataFrame()

    total_fips = data_all['fips_code']
    total_fips = np.unique(total_fips,return_counts=True)

    data_vis['FIPS'] = total_fips[0]
    data_vis['Total'] = total_fips[1]

    fl_school_fips_temp = np.empty(len(total_fips[0]),int)

    for f in range(0,len(total_fips[0])):
        if total_fips[0][f] in fl_school:
            fl_school_fips_temp[f] = fl_school[total_fips[0][f]]
        else:
            fl_school_fips_temp[f] = 0

    data_vis['Total school attendees'] = fl_school_fips_temp
    data_vis['Portion of school attendees'] = data_vis['Total school attendees'] / data_vis['Total']

    # EpiCast simulation analysis
    exposed_state = 1

    exposed_non_hispanic = []
    exposed_hispanic = []

    exposed_non_hispanic_overall_1 = []
    exposed_hispanic_overall_1 = []
    exposed_non_hispanic_overall_3 = []
    exposed_hispanic_overall_3 = []


    for e in events_files_all:
        if not SA:
            events_file = os.path.join(data_dir, "main_figures", f"run_{e:03}.events.bin")
        else:
            events_file = os.path.join(data_dir, "SA", f"run_{e:03}.events.bin")

        dn, dt, cn, events = read_event_file(events_file)
        if events_files_all[e] == 1 and not supplement: # check if model 1
            # Select school exposure only
            exposed = select_transitions(events, exposed_state, context_lo=2, context_hi=3) 

            full_ids = exposed['agent_id']
            all_states = agent_state(full_ids)

            # data = np.empty(len(full_ids), dtype=get_agent_type())

            ids = agent_id(full_ids[all_states == state])

            data = data_all.loc[ids]

            exposed_hispanic_temp = np.count_nonzero(data['person_ethnicity'] == 1)
            exposed_non_hispanic_temp = np.count_nonzero(data['person_ethnicity'] == 0)

            # Denominator: all hispanic/non-hispanic individuals
            # exposed_hispanic.append(exposed_hispanic_temp / hispanic_total)
            # exposed_non_hispanic.append(exposed_non_hispanic_temp / non_hispanic_total)

            # Denominator: all hispanic/non-hispanic school attendees
            exposed_hispanic.append(exposed_hispanic_temp / hispanic_school_total)
            exposed_non_hispanic.append(exposed_non_hispanic_temp / non_hispanic_school_total)

            # Measure school exposure in each census tract
            # exposed_fips = data['fips_code']
            # exposed_fips = np.unique(exposed_fips,return_counts=True)

            # exposed_dict = dict(zip(exposed_fips[0],exposed_fips[1]))

            # exposed_fips_tmp = np.empty(len(total_fips[0]), dtype=int)

            # for f in range(0,len(total_fips[0])):
            #     if total_fips[0][f] in exposed_fips[0]:
            #         exposed_fips_tmp[f] = exposed_dict[total_fips[0][f]]
            #     else:
            #         exposed_fips_tmp[f] = 0
        
            # data_vis['Exposed' + str(e)] = exposed_fips_tmp

            # Select all exposure events (regardless of context)
            exposed_overall = select_transitions(events, exposed_state, context_lo=-1, context_hi=255) 

            full_ids = exposed_overall['agent_id']
            all_states = agent_state(full_ids)

            ids = agent_id(full_ids[all_states == state])

            data = data_all.loc[ids]

            exposed_hispanic_temp = np.count_nonzero(data['person_ethnicity'] == 1)
            exposed_non_hispanic_temp = np.count_nonzero(data['person_ethnicity'] == 0)

            exposed_hispanic_overall_1.append(exposed_hispanic_temp / hispanic_total)
            exposed_non_hispanic_overall_1.append(exposed_non_hispanic_temp / non_hispanic_total)

            # Measure overall exposure in each census tract
            exposed_fips = data['fips_code']
            exposed_fips = np.unique(exposed_fips,return_counts=True)

            exposed_dict = dict(zip(exposed_fips[0],exposed_fips[1]))

            exposed_fips_tmp = np.empty(len(total_fips[0]), dtype=int)

            for f in range(0,len(total_fips[0])):
                if total_fips[0][f] in exposed_fips[0]:
                    exposed_fips_tmp[f] = exposed_dict[total_fips[0][f]]
                else:
                    exposed_fips_tmp[f] = 0
        
            data_vis['Exposed' + str(e)] = exposed_fips_tmp

        else: # model 3
            if supplement:
                # Select school exposure only
                exposed = select_transitions(events, exposed_state, context_lo=2, context_hi=3) 

                full_ids = exposed['agent_id']
                all_states = agent_state(full_ids)

                # data = np.empty(len(full_ids), dtype=get_agent_type())

                ids = agent_id(full_ids[all_states == state])

                data = data_all.loc[ids]

                exposed_hispanic_temp = np.count_nonzero(data['person_ethnicity'] == 1)
                exposed_non_hispanic_temp = np.count_nonzero(data['person_ethnicity'] == 0)

                # Denominator: all hispanic/non-hispanic individuals
                # exposed_hispanic.append(exposed_hispanic_temp / hispanic_total)
                # exposed_non_hispanic.append(exposed_non_hispanic_temp / non_hispanic_total)

                # Denominator: all hispanic/non-hispanic school attendees
                exposed_hispanic.append(exposed_hispanic_temp / hispanic_school_total)
                exposed_non_hispanic.append(exposed_non_hispanic_temp / non_hispanic_school_total)

                exposed_fips = data['fips_code']
                exposed_fips = np.unique(exposed_fips,return_counts=True)

                exposed_dict = dict(zip(exposed_fips[0],exposed_fips[1]))

                exposed_fips_tmp = np.empty(len(total_fips[0]), dtype=int)

                for f in range(0,len(total_fips[0])):
                    if total_fips[0][f] in exposed_fips[0]:
                        exposed_fips_tmp[f] = exposed_dict[total_fips[0][f]]
                    else:
                        exposed_fips_tmp[f] = 0
            
                data_vis['Exposed' + str(e)] = exposed_fips_tmp

            # Select all exposure events (regardless of context)
            exposed_overall = select_transitions(events, exposed_state, context_lo=-1, context_hi=255)

            full_ids = exposed_overall['agent_id']
            all_states = agent_state(full_ids)

            ids = agent_id(full_ids[all_states == state])

            data = data_all.loc[ids]

            exposed_hispanic_temp = np.count_nonzero(data['person_ethnicity'] == 1)
            exposed_non_hispanic_temp = np.count_nonzero(data['person_ethnicity'] == 0)

            exposed_hispanic_overall_3.append(exposed_hispanic_temp / hispanic_total)
            exposed_non_hispanic_overall_3.append(exposed_non_hispanic_temp / non_hispanic_total)

    # Aggregate to census tract level
    data_vis = data_vis.groupby(by=["FIPS"]).sum()

    filter_col = [col for col in data_vis if col.startswith('Exposed')]
    data_vis['Exposed'] = data_vis[filter_col].apply(np.mean,axis=1)

    # Denominator: all individuals in census tract
    data_vis['%Exposed'] = data_vis['Exposed'] / data_vis['Total']

    # Denominator: all school attendees in census tract
    # data_vis['%Exposed'] = data_vis['Exposed'] / data_vis['Total school attendees']

    plt.rcParams.update({'font.size': 16})

    if not supplement:
        fig, axes = plt.subplots(nrows=2, ncols=2)
        fig.set_figheight(11.25)
        fig.set_figwidth(15)
    else:
        fig, axes = plt.subplots(nrows=1, ncols=1)

    cmap = sb.color_palette("blend:#F8766D,#00BFC4", n_colors=2)

    if not supplement:
        # Figure 6A
        # fl_schoolgroups = pd.read_csv(os.path.join(data_dir, "aux", "fl_agent_schoolgroups.csv"))

        sg_size_dens_eth = pd.DataFrame()
        ethnicity_map = {0:'Non-Hisp.', 1:'Hisp.'}
        sg_size_dens_eth['Ethnicity'] = [ethnicity_map[h] for h in fl_schoolgroups['ethnicity']]
        sg_size_dens_eth['Schoolgroup size'] = fl_schoolgroups['sg_size']

        sb.kdeplot(data=sg_size_dens_eth, hue='Ethnicity', x='Schoolgroup size', palette=cmap, 
                common_norm=False, bw_adjust=2, clip=(1,45),ax=axes[0][0], fill=True, alpha=0.2)

        axes[0][0].text(
            0.0, 1.0, 'A)', transform=(
                axes[0][0].transAxes + ScaledTranslation(-20/72, +7/72, fig.dpi_scale_trans)),
            fontsize='medium', va='bottom', fontfamily='sans-serif')

        # Figure 6B

        fl_schoolgroups_fips = fl_schoolgroups.groupby(by='census_tracts').mean()

        data_vis.index = pd.to_numeric(data_vis.index)

        data_vis = data_vis.merge(fl_schoolgroups_fips,left_index=True,right_index=True)

        data_vis['Mean schoolgroup size'] = data_vis['sg_size']

        # cmap_vrid = sb.color_palette("viridis")

        sb.scatterplot(data=data_vis, x="Portion of school attendees", y="%Exposed", hue="Mean schoolgroup size", palette='viridis', alpha=0.5,  ax=axes[0][1])
        # axes[0][1].set_xlabel('Mean schoolgroup size among students/teachers in tract')
        axes[0][1].set_xlabel('Proportion living in tract attending school')

        # Portion of individuals infected at school in census tract
        # axes[0][1].set_ylabel('Proportion infected w/in school')

        # Portion of individuals infected
        axes[0][1].set_ylabel('Proportion infected overall')

        # Portion of school attendees infected at school in census tract
        # axes[0][1].set_ylabel('Proportion of school attendees infected \nw/in school')

        axes[0][1].text(
            0.0, 1.0, 'B)', transform=(
                axes[0][1].transAxes + ScaledTranslation(-20/72, +7/72, fig.dpi_scale_trans)),
            fontsize='medium', va='bottom', fontfamily='sans-serif')


    # Figure 6C
    six_c_axes = axes if supplement else axes[1][0]

    school_exp_groups = pd.DataFrame()
    school_exp_groups['Racial/ethnic group'] = ['Non-Hisp.']*len(exposed_non_hispanic) + ['Hisp.']*len(exposed_non_hispanic)
    school_exp_groups['Exposure'] = exposed_non_hispanic + exposed_hispanic

    # print(school_exp_groups.groupby(by=['Racial/ethnic group']).mean())
    # print(school_exp_groups.groupby(by=['Racial/ethnic group']).var())

    sb.barplot(school_exp_groups,x='Racial/ethnic group',
               y='Exposure',errorbar=None, palette=cmap,ax=six_c_axes)
    sb.stripplot(school_exp_groups,x='Racial/ethnic group',
                 y='Exposure', alpha=0.4, s=10, color='black',ax=six_c_axes)
    
    # Denominator: all individuals in racial/ethnic group
    # six_c_axes.set_ylabel('Proportion infected w/in school')

    # Denominator: all school attendees in racial/ethnic group
    six_c_axes.set_ylabel('Proportion of school attendees infected \nw/in school')

    if not supplement:
        six_c_axes.text(
            0.0, 1.0, 'C)', transform=(
                six_c_axes.transAxes + ScaledTranslation(-20/72, +7/72, fig.dpi_scale_trans)),
            fontsize='medium', va='bottom', fontfamily='sans-serif')

    # Figure 6D
    if not supplement:

        # Overall infections under each model
        school_exp_groups = pd.DataFrame()
        school_exp_groups['Racial/ethnic group'] = ['Non-Hisp.']*len(exposed_non_hispanic_overall_1) + ['Hisp.']*len(exposed_non_hispanic_overall_1) + ['Non-Hisp.']*len(exposed_non_hispanic_overall_3) + ['Hisp.']*len(exposed_non_hispanic_overall_3)
        school_exp_groups['Model'] = ['Model 1']*len(exposed_non_hispanic_overall_1)*2 + ['Model 3']*len(exposed_non_hispanic_overall_1)*2
        school_exp_groups['Exposure'] = exposed_non_hispanic_overall_1 + exposed_hispanic_overall_1 + exposed_non_hispanic_overall_3 + exposed_hispanic_overall_3

        # print(school_exp_groups.groupby(by=['Racial/ethnic group','Model']).mean())
        # print(school_exp_groups.groupby(by=['Racial/ethnic group','Model']).var())

        # sb.barplot(school_exp_groups,x='Racial/ethnic group',y='Exposure', hue='Model', errorbar=None, palette="YlGnBu", ax=axes[1][1])
        # sb.stripplot(school_exp_groups,x='Racial/ethnic group',y='Exposure', hue='Model', s=10, color='black', dodge=True, alpha=0.4, ax=axes[1][1])

        # axes[1][1].set_xlabel('Racial/ethnic group')
        # axes[1][1].set_ylabel('Proportion infected overall')

        # Overall difference in infection rate between hispanic & non-hispanic pop.
        school_exp_groups = pd.DataFrame()
        school_exp_groups['Model'] = ['Model 1']*len(exposed_non_hispanic_overall_1) + ['Model 3']*len(exposed_non_hispanic_overall_1)
        school_exp_groups['Exposure'] = [ exposed_hispanic_overall_1[h] - exposed_non_hispanic_overall_1[h] for h in range(len(exposed_non_hispanic_overall_1))] + [ exposed_hispanic_overall_3[h] - exposed_non_hispanic_overall_3[h] for h in range(len(exposed_non_hispanic_overall_3))]

        sb.barplot(school_exp_groups,x='Model',
                y='Exposure', errorbar=None, 
                palette="YlGnBu", ax=axes[1][1])
        sb.stripplot(school_exp_groups,x='Model',
                    y='Exposure', alpha=0.5, s=10, 
                    color='black', ax=axes[1][1], dodge=True)

        axes[1][1].set_xlabel('Model')
        # axes[1][1].set_ylabel('Overall difference (%Hispanic pop. infected - \n%non-Hispanic pop. infected)')
        axes[1][1].set_ylabel('Δ proportion infected overall\n(Hispanic - non-Hispanic)')


        # handles, labels = axes[1][1].get_legend_handles_labels()
        # axes[1][1].legend(handles[0:2], labels[0:2])

        axes[1][1].text(
            0.0, 1.0, 'D)', transform=(
                axes[1][1].transAxes + ScaledTranslation(-20/72, +7/72, fig.dpi_scale_trans)),
            fontsize='medium', va='bottom', fontfamily='sans-serif')


    six_c_axes.spines[['right','top']].set_visible(False)
    if not supplement:
        axes[0][0].spines[['right','top']].set_visible(False)
        axes[0][1].spines[['right','top']].set_visible(False)
        # axes[1][0].spines[['right','top']].set_visible(False)
        axes[1][1].spines[['right','top']].set_visible(False)


    plt.tight_layout()
    if supplement:
        plt.savefig('figS2.pdf')
    elif SA:
        plt.savefig('figS6.pdf')
    else:
        plt.savefig('fig6.pdf')

    return
# ============================================================================ #
def figure_s1_gen(data_dir):
    figure_4_gen(data_dir, 35, "ethnicity")
    os.rename("fig4_35_ethnicity.pdf", "figS1.pdf")
# ============================================================================ #
def figure_s2_gen(data_dir):
    figure_6_gen(data_dir, supplement=True)
# ============================================================================ #
def figure_s3_gen(data_dir):
    figure_4_gen(data_dir, 12, "ethnicity")
    os.rename("fig4_12_ethnicity.pdf", "figS3.pdf")
# ============================================================================ #
def figure_s4_gen(data_dir):
    figure_4_gen(data_dir, SA=True)
# ============================================================================ #
def figure_s5_gen(data_dir):
    figure_5_gen(data_dir, SA=True)
# ============================================================================ #
def figure_s6_gen(data_dir):
    figure_6_gen(data_dir, SA=True)
# ============================================================================ #
if __name__ == "__main__":

    if len(sys.argv) < 3:
        print("\n[ERROR]: Not enough inputs, at least 2 are required - see below:\n")
        print("Usage: epicast_heterogeneity.py <fig_num> <data_dir>")
        print("   <fig_num>: one of [4,5,6,s1,s2,s3]")
        print("   <data_dir>: full path to directory containing data for the paper")
        sys.exit(1)

    fig_num = sys.argv[1]
    data_dir = sys.argv[2]

    if fig_num.lower() == "s1":
        figure_s1_gen(data_dir)

    elif fig_num.lower() == "s2":
        figure_s2_gen(data_dir)

    elif fig_num.lower() == "s3":
        figure_s3_gen(data_dir)

    elif fig_num.lower() == "s4":
        figure_s4_gen(data_dir)

    elif fig_num.lower() == "s5":
        figure_s5_gen(data_dir)

    elif fig_num.lower() == "s6":
        figure_s6_gen(data_dir)

    elif fig_num == "4":
        # default to NM/race if no state/variable is specified
        state = 35 if len(sys.argv) < 4 else int(sys.argv[3])
        compare = "race" if len(sys.argv) < 5 else sys.argv[4]

        figure_4_gen(data_dir, state, compare)

        # regular figure 4
        if state == 35 and compare == "race":
            os.rename("fig4_35_race.pdf", "fig4.pdf")

    elif fig_num == "5":
        figure_5_gen(data_dir)

    elif fig_num == "6":
        figure_6_gen(data_dir)

    else:
        raise Exception(f"Figure number {fig_num} is not valid, must be in [4,5,6,s1,s2,s3]")