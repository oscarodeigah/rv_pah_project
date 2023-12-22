import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib as mpl

import numpy as np
import json
import yaml

#sns.set_theme()

def plot_data_assimilation_results(pv_results, active_results):

    mpl.rcParams['lines.linewidth'] = 2

    with open(pv_results) as json_file:
        pv_res = json.load(json_file)

    with open(active_results) as json_file:
        active_res = json.load(json_file)

    fig, axs = plt.subplots(nrows=2, ncols=4, figsize=(20,14))
    rc('font', weight='bold')

    axs[0,0].plot(pv_res['CNT']['simulated_volume_rv'][1:], pv_res['CNT']['rv_pressure'][1:], 'k-')
    axs[0,0].plot(pv_res['CNT']['measured_volume_rv'][1:], pv_res['CNT']['rv_pressure'][1:], 'ro')
    axs[0,0].set_title('Control', fontsize=20, fontweight='bold')

    axs[0,1].plot(pv_res['PAHw4']['simulated_volume_rv'][1:], pv_res['PAHw4']['rv_pressure'][1:], 'k-')
    axs[0,1].plot(pv_res['PAHw4']['measured_volume_rv'][1:], pv_res['PAHw4']['rv_pressure'][1:], 'ro')
    axs[0,1].set_title('SuHx Week 4', fontsize=20, fontweight='bold')

    axs[0,2].plot(pv_res['PAHw8']['simulated_volume_rv'][1:], pv_res['PAHw8']['rv_pressure'][1:], 'k-')
    axs[0,2].plot(pv_res['PAHw8']['measured_volume_rv'][1:], pv_res['PAHw8']['rv_pressure'][1:], 'ro')
    axs[0,2].set_title('SuHx Week 8', fontsize=20, fontweight='bold')

    axs[0,3].plot(pv_res['PAHw12']['simulated_volume_rv'][1:], pv_res['PAHw12']['rv_pressure'][1:], 'k-')
    axs[0,3].plot(pv_res['PAHw12']['measured_volume_rv'][1:], pv_res['PAHw12']['rv_pressure'][1:], 'ro')
    axs[0,3].set_title('SuHx Week 12', fontsize=20, fontweight='bold')

    axs[1,0].plot(active_res['CNT'])
    axs[1,1].plot(active_res['PAHw4'])
    axs[1,2].plot(active_res['PAHw8'])
    axs[1,3].plot(active_res['PAHw12'])

    for ax in axs.flat:
        #ax.grid(True)
        ax.tick_params(axis='both', which='major', labelsize=18)
        ax.tick_params(axis='both', which='minor', labelsize=18)
        ax.spines['bottom'].set_linewidth(2)
        ax.spines['top'].set_linewidth(2)
        ax.spines['left'].set_linewidth(2)
        ax.spines['right'].set_linewidth(2)

    # Add a legend outside and below the plot
    #fig.legend(['Simulated', 'Measured'], loc='center right', ncol=1, fontsize=14)
    fig.subplots_adjust(hspace=0.3, bottom=0.1)

    # Draw vertical lines through the ED and ES points of the bottom row
    x_points = [3, 23]  # coordinates of the ED and ES points respectively
    for ax in (axs[1,0], axs[1,1], axs[1,2], axs[1,3]):
        for p in x_points:
            ax.axvline(p, color='black', linestyle='dotted', linewidth=2)

    # Add a general x label to the top row
    fig.text(0.51, 0.49, 'Volume (uL)', ha='center', va='center', fontsize=20)

    # Add a general y label to the top row
    fig.text(0.09, 0.7, 'Pressure (kPa)', ha='center', va='center', rotation='vertical', fontsize=20)

    # Add a general y label to the bottom row
    fig.text(0.09, 0.27, 'Active stress, $T_{a}$ (kPa)', ha='center', va='center', rotation='vertical', fontsize=20)

    # Hide the y labels and yticks for the second to last columns
    #for ax in axs.flat:
    #    ax.label_outer()

    for ax in (axs[0,0], axs[0,1], axs[0,2]):
        ax.set_ylim(0, 9)
        ax.set_xlim(50, 250)
        ax.set_yticks([0, 2, 4, 6, 8])
        ax.set_xticklabels([50, "", 150, "", 250], fontweight='bold')

    axs[0,0].set_yticklabels([0, 2, 4, 6, 8], fontweight='bold')
    axs[0,3].set_xlim(150, 350)
    axs[0,3].set_ylim(0, 9)
    axs[0,3].set_yticks([0, 2, 4, 6, 8])
    axs[0,3].set_yticklabels([0, 5, 10, 15, 20], fontweight='bold')
    axs[0,3].set_xticklabels([150, "", 250, "", 350], fontweight='bold')

    for ax in (axs[0,1], axs[0,2], axs[0,3]):
        ax.set_yticklabels(["", "", "", "", ""])

    for ax in (axs[1,0], axs[1,1], axs[1,2], axs[1,3]):
        ax.set_ylim(-0.5, 40.0)
        ax.set_yticks([0, 10, 20, 30, 40])
        ax.set_xticks([3, 23])

    axs[1,0].set_yticklabels([0, 10, 20, 30, 40], fontweight='bold')

    for ax in (axs[1,1], axs[1,2], axs[1,3]):
        ax.set_yticklabels(["", "", "", "", ""])
    for ax in (axs[1,0], axs[1,1], axs[1,2], axs[1,3]):
        a = ax.get_xticks().tolist()
        a[0] = 'ED'; a[1] = 'ES'; ax.set_xticklabels(a, fontweight='bold')

    fig.savefig('data_assimilation_results.pdf', dpi=300)


def plot_mesh_convergence_analysis_stress_strain_results(stress_results, strain_results):

    mpl.rcParams['lines.linewidth'] = 3

    with open(stress_results) as json_file:
        stress_res = json.load(json_file)

    with open(strain_results) as json_file:
        strain_res = json.load(json_file)

    fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(20,14))
    rc('font', weight='bold')

    axs[0,0].plot(stress_res['lowRes']['fiber'])
    axs[0,0].plot(stress_res['midlowRes']['fiber'])
    axs[0,0].plot(stress_res['midhighRes']['fiber'])
    axs[0,0].plot(stress_res['highRes']['fiber'])
    axs[0,0].set_ylabel('Fiber stress (kPa)', fontsize=17, fontweight='bold')
    fig.text(0.24, 0.91, 'Fiber direction', ha='center', va='center', fontsize=19, fontweight='bold')

    axs[0,1].plot(stress_res['lowRes']['circumferential'])
    axs[0,1].plot(stress_res['midlowRes']['circumferential'])
    axs[0,1].plot(stress_res['midhighRes']['circumferential'])
    axs[0,1].plot(stress_res['highRes']['circumferential'])
    axs[0,1].set_ylabel('Circumferential stress (kPa)', fontsize=17, fontweight='bold')
    fig.text(0.51, 0.91, 'Circumferential direction', ha='center', va='center', fontsize=19, fontweight='bold')

    axs[0,2].plot(stress_res['lowRes']['longitudinal'])
    axs[0,2].plot(stress_res['midlowRes']['longitudinal'])
    axs[0,2].plot(stress_res['midhighRes']['longitudinal'])
    axs[0,2].plot(stress_res['highRes']['longitudinal'])
    axs[0,2].set_ylabel('AOT stress (kPa)', fontsize=17, fontweight='bold')
    fig.text(0.79, 0.91, 'AOT direction', ha='center', va='center', fontsize=19, fontweight='bold')

    axs[1,0].plot(strain_res['lowRes']['fiber'])
    axs[1,0].plot(strain_res['midlowRes']['fiber'])
    axs[1,0].plot(strain_res['midhighRes']['fiber'])
    axs[1,0].plot(strain_res['highRes']['fiber'])
    axs[1,0].set_ylabel('Fiber strain', fontsize=17, fontweight='bold')

    axs[1,1].plot(strain_res['lowRes']['circumferential'])
    axs[1,1].plot(strain_res['midlowRes']['circumferential'])
    axs[1,1].plot(strain_res['midhighRes']['circumferential'])
    axs[1,1].plot(strain_res['highRes']['circumferential'])
    axs[1,1].set_ylabel('Circumferential strain', fontsize=17, fontweight='bold')

    axs[1,2].plot(strain_res['lowRes']['longitudinal'])
    axs[1,2].plot(strain_res['midlowRes']['longitudinal'])
    axs[1,2].plot(strain_res['midhighRes']['longitudinal'])
    axs[1,2].plot(strain_res['highRes']['longitudinal'])
    axs[1,2].set_ylabel('AOT strain', fontsize=17, fontweight='bold')

    # Add a legend outside and below the plot
    fig.legend(['Low resolution mesh', 'Medium-low resolution mesh', 'Medium-high resolution mesh', 'High resolution mesh'], loc='lower center', ncol=4, fontsize=17)
    fig.subplots_adjust(hspace=0.15, bottom=0.1)

    # Draw vertical lines through the ED and ES points
    x_points = [3, 23]  # coordinates of the ED and ES points respectively
    for ax in axs.flat:
        for p in x_points:
            ax.axvline(p, color='black', linestyle='dotted', linewidth=2)

    # add grid lines to all axes
    for ax in axs.flat:
        #ax.grid(True)
        ax.tick_params(axis='both', which='major', labelsize=16)
        ax.tick_params(axis='both', which='minor', labelsize=16)

    for ax in (axs[0,0], axs[0,1], axs[0,2]):
        ax.set_ylim(0.0, 22.0)
        ax.set_yticks([0, 5, 10, 15, 20])
        ax.set_xticks([3, 23])

    axs[0,0].set_yticklabels([0, 5, 10, 15, 20], fontweight='bold')

    for ax in (axs[0,1], axs[0,2]):
        ax.set_yticklabels(["", "", "", "", ""])

    for ax in (axs[0,0], axs[0,1], axs[0,2]):
        ax.set_xticklabels(["", ""])

    for ax in (axs[1,0], axs[1,1], axs[1,2]):
        ax.set_ylim(-0.05, 0.55)
        ax.set_yticks([0, 0.1, 0.2, 0.3, 0.4, 0.5])
        ax.set_xticks([3, 23])

    axs[1,0].set_yticklabels([0, 0.1, 0.2, 0.3, 0.4, 0.5], fontweight='bold')

    for ax in (axs[1,1], axs[1,2]):
        ax.set_yticklabels(["", "", "", "", "", ""])

    #for ax in axs.flat:
    #    a = ax.get_xticks().tolist()
    #    a[0] = 'ED'; a[1] = 'ES'; ax.set_xticklabels(a, fontweight='bold')

    for ax in (axs[1,0], axs[1,1], axs[1,2]):
        a = ax.get_xticks().tolist()
        a[0] = 'ED'; a[1] = 'ES'; ax.set_xticklabels(a, fontweight='bold')

    for ax in axs.flat:
        ax.spines['bottom'].set_linewidth(2)
        ax.spines['top'].set_linewidth(2)
        ax.spines['left'].set_linewidth(2)
        ax.spines['right'].set_linewidth(2)

    # Save the figure
    fig.savefig('mesh_convergence_analysis_stress_strain_results.pdf', dpi=300)


def plot_stress_strain_results(stress_results, strain_results):

    mpl.rcParams['lines.linewidth'] = 3

    with open(stress_results) as json_file:
        stress_res = json.load(json_file)

    with open(strain_results) as json_file:
        strain_res = json.load(json_file)

    fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(20,14))
    rc('font', weight='bold')

    axs[0,0].plot(stress_res['CNT']['fiber'])
    axs[0,0].plot(stress_res['PAHw4']['fiber'])
    axs[0,0].plot(stress_res['PAHw8']['fiber'])
    axs[0,0].plot(stress_res['PAHw12']['fiber'])
    axs[0,0].set_ylabel('Fiber stress (kPa)', fontsize=17, fontweight='bold')
    #axs[0,0].set_title('Fiber direction', fontsize=17, fontweight='bold')
    fig.text(0.24, 0.91, 'Fiber direction', ha='center', va='center', fontsize=19, fontweight='bold')

    axs[0,1].plot(stress_res['CNT']['circumferential'])
    axs[0,1].plot(stress_res['PAHw4']['circumferential'])
    axs[0,1].plot(stress_res['PAHw8']['circumferential'])
    axs[0,1].plot(stress_res['PAHw12']['circumferential'])
    axs[0,1].set_ylabel('Circumferential stress (kPa)', fontsize=17, fontweight='bold')
    #axs[0,1].set_title('Circumferential direction', fontsize=17, fontweight='bold')
    fig.text(0.51, 0.91, 'Circumferential direction', ha='center', va='center', fontsize=19, fontweight='bold')

    axs[0,2].plot(stress_res['CNT']['longitudinal'])
    axs[0,2].plot(stress_res['PAHw4']['longitudinal'])
    axs[0,2].plot(stress_res['PAHw8']['longitudinal'])
    axs[0,2].plot(stress_res['PAHw12']['longitudinal'])
    axs[0,2].set_ylabel('AOT stress (kPa)', fontsize=17, fontweight='bold')
    #axs[0,2].set_title('AOT direction', fontsize=17, fontweight='bold')
    fig.text(0.79, 0.91, 'AOT direction', ha='center', va='center', fontsize=19, fontweight='bold')

    axs[1,0].plot(strain_res['CNT']['fiber'])
    axs[1,0].plot(strain_res['PAHw4']['fiber'])
    axs[1,0].plot(strain_res['PAHw8']['fiber'])
    axs[1,0].plot(strain_res['PAHw12']['fiber'])
    axs[1,0].set_ylabel('Fiber strain', fontsize=17, fontweight='bold')

    axs[1,1].plot(strain_res['CNT']['circumferential'])
    axs[1,1].plot(strain_res['PAHw4']['circumferential'])
    axs[1,1].plot(strain_res['PAHw8']['circumferential'])
    axs[1,1].plot(strain_res['PAHw12']['circumferential'])
    axs[1,1].set_ylabel('Circumferential strain', fontsize=17, fontweight='bold')

    axs[1,2].plot(strain_res['CNT']['longitudinal'])
    axs[1,2].plot(strain_res['PAHw4']['longitudinal'])
    axs[1,2].plot(strain_res['PAHw8']['longitudinal'])
    axs[1,2].plot(strain_res['PAHw12']['longitudinal'])
    axs[1,2].set_ylabel('AOT strain', fontsize=17, fontweight='bold')

    # Add a legend outside and below the plot
    fig.legend(['Control', 'SuHx Week 4', 'SuHx Week 8', 'SuHx Week 12'], loc='lower center', ncol=4, fontsize=17)
    fig.subplots_adjust(hspace=0.15, bottom=0.1)

    # Draw vertical lines through the ED and ES points
    x_points = [3, 23]  # coordinates of the ED and ES points respectively
    for ax in axs.flat:
        for p in x_points:
            ax.axvline(p, color='black', linestyle='dotted', linewidth=2)

    # add grid lines to all axes
    for ax in axs.flat:
        #ax.grid(True)
        ax.tick_params(axis='both', which='major', labelsize=16)
        ax.tick_params(axis='both', which='minor', labelsize=16)
        #ax.spines['top'].set_visible(False)
        #ax.spines['right'].set_visible(False)

    for ax in (axs[0,0], axs[0,1], axs[0,2]):
        ax.set_ylim(-0.5, 23.0)
        ax.set_yticks([0, 5, 10, 15, 20])
        ax.set_xticks([3, 23])

    axs[0,0].set_yticklabels([0, 5, 10, 15, 20], fontweight='bold')

    for ax in (axs[0,1], axs[0,2]):
        ax.set_yticklabels(["", "", "", "", ""])

    for ax in (axs[0,0], axs[0,1], axs[0,2]):
        ax.set_xticklabels(["", ""])

    for ax in (axs[1,0], axs[1,1], axs[1,2]):
        ax.set_ylim(-0.03, 0.52)
        ax.set_yticks([0, 0.1, 0.2, 0.3, 0.4, 0.5])
        ax.set_xticks([3, 23])

    axs[1,0].set_yticklabels([0, 0.1, 0.2, 0.3, 0.4, 0.5], fontweight='bold')

    for ax in (axs[1,1], axs[1,2]):
        ax.set_yticklabels(["", "", "", "", "", ""])

    #for ax in axs.flat:
    #    a = ax.get_xticks().tolist()
    #    a[0] = 'ED'; a[1] = 'ES'; ax.set_xticklabels(a, fontweight='bold')

    for ax in (axs[1,0], axs[1,1], axs[1,2]):
        a = ax.get_xticks().tolist()
        a[0] = 'ED'; a[1] = 'ES'; ax.set_xticklabels(a, fontweight='bold')

    for ax in axs.flat:
        ax.spines['bottom'].set_linewidth(2)
        ax.spines['top'].set_linewidth(2)
        ax.spines['left'].set_linewidth(2)
        ax.spines['right'].set_linewidth(2)

    # Save the figure
    fig.savefig('stress_strain_results.pdf', dpi=300)


def plot_effect_of_remodeling():

    fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(24,16))
    rc('font', weight='bold')

    for ax in axs.flat:
        ax.tick_params(axis='both', which='major', labelsize=18)
        ax.tick_params(axis='both', which='minor', labelsize=18)
        ax.spines['bottom'].set_linewidth(2.5)
        ax.spines['top'].set_linewidth(2.5)
        ax.spines['left'].set_linewidth(2.5)
        ax.spines['right'].set_linewidth(2.5)

    width = 0.2
    x = np.arange(3)

    # ED fiber stress [w4, w8, w12]
    y00_1 = [2.42, 4.76, 8.47] # No remodeling 
    y00_2 = [2.18, 4.12, 7.48] # No geometric remodeling
    y00_3 = [1.09, 1.85, 3.03] # No material remodeling
    y00_4 = [0.95, 1.58, 2.7] # Geometric & material remodeling

    # ES fiber stress [w4, w8, w12]
    y01_1 = [52.6, 85.81, 74.11] # No remodeling 
    y01_2 = [49.61, 75.19, 71.53] # No geometric remodeling
    y01_3 = [16.46, 20.37, 18.57] # No material remodeling
    y01_4 = [16.38, 20.33, 20.78] # Geometric & material remodelin18
    
    # ED fiber strain [w4, w8, w12]
    y10_1 = [0.27, 0.34, 0.39] # No remodeling 
    y10_2 = [0.22, 0.27, 0.34] # No geometric remodeling
    y10_3 = [0.17, 0.21, 0.23] # No material remodeling
    y10_4 = [0.14, 0.16, 0.2] # Geometric & material remodeling

    # ES fiber strain [w4, w8, w12]
    y11_1 = [0.31, 0.45, 0.41] # No remodeling 
    y11_2 = [0.29, 0.38, 0.45] # No geometric remodeling
    y11_3 = [-0.001, 0.024, -0.01] # No material remodeling
    y11_4 = [0.01, 0.02, 0.1] # Geometric & material remodeling


    axs[0,0].bar(x-0.2, y00_1, width, edgecolor='black', color='lavenderblush')
    axs[0,0].bar(x, y00_2, width, edgecolor='black', color='palevioletred')
    axs[0,0].bar(x+0.2, y00_3, width, edgecolor='black', color='crimson')
    axs[0,0].bar(x+0.4, y00_4, width, edgecolor='black', color='black')#, hatch='\/')
    #axs[0,0].set_ylabel('ED Fiber Stress (kPa)', fontsize=18, fontweight='bold')
    axs[0,0].axhline(0.36, color="dimgray", linewidth=3)
    axs[0,0].text(1.005, 0.36, "0.36", va='center', ha="left", bbox=dict(facecolor="w",alpha=0.5), transform=axs[0,0].get_yaxis_transform(), fontsize=18)

    axs[0,1].bar(x-0.2, y01_1, width, edgecolor='black', color='lavenderblush')
    axs[0,1].bar(x, y01_2, width, edgecolor='black', color='palevioletred')
    axs[0,1].bar(x+0.2, y01_3, width, edgecolor='black', color='crimson')
    axs[0,1].bar(x+0.4, y01_4, width, edgecolor='black', color='black')#, hatch='\/')
    #axs[0,1].set_ylabel('ES Fiber Stress (kPa)', fontsize=18, fontweight='bold')
    axs[0,1].axhline(14.48, color="dimgray", linewidth=3)
    axs[0,1].text(1.005, 14.48, "14.48", va='center', ha="left", bbox=dict(facecolor="w",alpha=0.5), transform=axs[0,1].get_yaxis_transform(), fontsize=18)

    axs[1,0].bar(x-0.2, y10_1, width, edgecolor='black', color='lavenderblush')
    axs[1,0].bar(x, y10_2, width, edgecolor='black', color='palevioletred')
    axs[1,0].bar(x+0.2, y10_3, width, edgecolor='black', color='crimson')
    axs[1,0].bar(x+0.4, y10_4, width, edgecolor='black', color='black')#, hatch='\/')
    #axs[1,0].set_ylabel('ED Fiber Strain', fontsize=18, fontweight='bold')
    axs[1,0].axhline(0.14, color="dimgray", linewidth=3)
    axs[1,0].text(1.005, 0.14, "0.14", va='center', ha="left", bbox=dict(facecolor="w",alpha=0.5), transform=axs[1,0].get_yaxis_transform(), fontsize=18)

    axs[1,1].bar(x-0.2, y11_1, width, edgecolor='black', color='lavenderblush')
    axs[1,1].bar(x, y11_2, width, edgecolor='black', color='palevioletred')
    axs[1,1].bar(x+0.2, y11_3, width, edgecolor='black', color='crimson')
    axs[1,1].bar(x+0.4, y11_4, width, edgecolor='black', color='black')#, hatch='\/')
    #axs[1,1].set_ylabel('ES Fiber Strain', fontsize=18, fontweight='bold')
    axs[1,1].axhline(0.0002, color="dimgray", linewidth=3)
    axs[1,1].text(1.005, 0.0002, "0.0002", va='center', ha="left", bbox=dict(facecolor="w",alpha=0.5), transform=axs[1,1].get_yaxis_transform(), fontsize=18)

    fig.legend(['Control value', 'No remodeling', 'Only material remodeling', 'Only geometric remodeling', 'Both geometric & material remodeling'], loc='lower center', ncol=5, fontsize=16)
    fig.subplots_adjust(hspace=0.25, wspace=0.3, bottom=0.1)

    # Add label to axs[0,0]
    fig.text(0.3, 0.49, 'Week', ha='center', va='center', fontsize=18)
    fig.text(0.09, 0.7, 'ED fiber stress (kPa)', ha='center', va='center', rotation='vertical', fontsize=18)

    # Add label to axs[0,1]
    fig.text(0.74, 0.49, 'Week', ha='center', va='center', fontsize=18)
    fig.text(0.527, 0.7, 'ES fiber stress (kPa)', ha='center', va='center', rotation='vertical', fontsize=18)

    # Add label to axs[1,0]
    fig.text(0.3, 0.06, 'Week', ha='center', va='center', fontsize=18)
    fig.text(0.09, 0.27, 'ED fiber strain', ha='center', va='center', rotation='vertical', fontsize=18)

    # Add label to axs[1,1]
    fig.text(0.74, 0.06, 'Week', ha='center', va='center', fontsize=18)
    fig.text(0.527, 0.27, 'ES fiber strain', ha='center', va='center', rotation='vertical', fontsize=18)

    # axs[0,0]
    axs[0,0].set_ylim(0, 10)
    axs[0,0].set_yticks([0, 2, 4, 6, 8, 10])
    axs[0,0].set_yticklabels([0, 2, 4, 6, 8, 10], fontweight='bold')

    # axs[0,1]
    axs[0,1].set_ylim(0, 100)
    axs[0,1].set_yticks([0, 20, 40, 60, 80, 100])
    axs[0,1].set_yticklabels([0, 20, 40, 60, 80, 100], fontweight='bold')

    # axs[1,0]
    axs[1,0].set_ylim(0, 0.5)
    axs[1,0].set_yticks([0, 0.1, 0.2, 0.3, 0.4, 0.5])
    axs[1,0].set_yticklabels([0, 0.1, 0.2, 0.3, 0.4, 0.5], fontweight='bold')

    # axs[1,1]
    axs[1,1].set_ylim(-0.02, 0.5)
    axs[1,1].set_yticks([0, 0.1, 0.2, 0.3, 0.4, 0.5])
    axs[1,1].set_yticklabels([0, 0.1, 0.2, 0.3, 0.4, 0.5], fontweight='bold')

    for ax in axs.flat:
        ax.set_xticks(x, ['4', '8', '12'])
        ax.set_xticklabels(['4', '8', '12'], fontweight='bold')
        #ax.set_xlabel('Week', fontsize=12)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    # Save the figure
    fig.savefig('effect_of_geo_remodeling_new_test.pdf', dpi=300)


def plot_lv_data(lv_pv_data, lv_active_data):

    mpl.rcParams['lines.linewidth'] = 3

    fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(14,6))
    rc('font', weight='bold')

    with open(lv_pv_data) as f:
        pv_data = yaml.safe_load(f)
    
    LVP = pv_data['LVP']; LVV = pv_data['LVV']

    with open(lv_active_data) as f:
        lv_active = json.load(f)

    axs[0].plot(LVV, LVP, 'ro')
 
    axs[0].tick_params(axis='both', which='major', labelsize=18)
    axs[0].tick_params(axis='both', which='minor', labelsize=18)
    axs[0].spines['bottom'].set_linewidth(2)
    axs[0].spines['top'].set_linewidth(2)
    axs[0].spines['left'].set_linewidth(2)
    axs[0].spines['right'].set_linewidth(2)
    axs[0].set_ylim(0,15)
    axs[0].set_yticks([0,2.5,5,7.5,10,12.5,15])
    axs[0].set_yticklabels(["0", "", "5", "", "10", "", "15"], fontweight='bold')
    axs[0].set_ylabel('Pressure (kPa)', fontweight='bold', fontsize=20)
    axs[0].set_xlim(60,240)
    axs[0].set_xticks([60,90,120,150,180,210,240])
    axs[0].set_xticklabels(["60", "", "120", "", "180", "", "240"], fontweight='bold')
    axs[0].set_xlabel('Volume (uL)', fontweight='bold', fontsize=20)

    axs[1].plot(lv_active["LV_active_data"])

    x_points = [3, 23]
    for p in x_points:
            axs[1].axvline(p, color='black', linestyle='dotted', linewidth=2)

    axs[1].set_xticks([3, 23])
    a = axs[1].get_xticks().tolist()
    a[0] = 'ED'; a[1] = 'ES'; axs[1].set_xticklabels(a, fontweight='bold')

    axs[1].set_ylim(-0.5,110)
    axs[1].set_yticks([0,20,40,60,80,100])
    axs[1].set_yticklabels(["0", "20", "40", "60", "80", "100"], fontweight='bold')
    axs[1].set_ylabel('Active stress, $\mathbf{T_{a}}$ (kPa)', fontweight='bold', fontsize=20)

    axs[1].tick_params(axis='both', which='major', labelsize=18)
    axs[1].tick_params(axis='both', which='minor', labelsize=18)
    axs[1].spines['bottom'].set_linewidth(2)
    axs[1].spines['top'].set_linewidth(2)
    axs[1].spines['left'].set_linewidth(2)
    axs[1].spines['right'].set_linewidth(2)

    fig.subplots_adjust(wspace=0.3)

    fig.savefig('lv_data.pdf', dpi=300)


if __name__ == '__main__':
    plot_data_assimilation_results('all_PV_results.json', 'all_active_results.json')
    plot_lv_data("../data/sample_datafiles/CNT.yml", "../data/sample_datafiles/LV_active_data.json")
    plot_mesh_convergence_analysis_stress_strain_results('MeshConvergenceAnalysis_all_stress_results.json', 'MeshConvergenceAnalysis_all_strain_results.json')
    plot_stress_strain_results('all_stress_results.json', 'all_strain_results.json')
    plot_effect_of_remodeling()
    
    
    
