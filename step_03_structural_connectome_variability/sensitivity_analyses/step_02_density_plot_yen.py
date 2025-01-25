import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

# Set global font settings at the start
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 36
plt.rcParams['axes.labelcolor'] = 'black'
plt.rcParams['xtick.color'] = 'black'
plt.rcParams['ytick.color'] = 'black'

def create_cmp():
    from matplotlib import cm
    from matplotlib.colors import ListedColormap
    blue = cm.get_cmap('Blues_r', 256)
    newcolors = blue(np.linspace(0.5, 0.85, 256))
    newcmp = ListedColormap(newcolors)
    return newcmp

def density_show(input_file,x_label,y_label,out_path):
    cmap = create_cmp()
    data = pd.read_csv(input_file)
    x, y = 'mat_a','mat_b'
    fig, ax = plt.subplots(figsize=(11, 12))
    
    # Plot with consistent font settings
    sns.kdeplot(data=data, x=x, y=y, fill=True, levels=12, thresh=0.1, cmap=cmap, ax=ax)
    sns.regplot(data=data, x=x, y=y, scatter=False, color='steelblue', 
                line_kws={'linewidth': '8'})
    
    # Axis and spine settings
    ax.spines.right.set_visible(False)
    ax.spines.top.set_visible(False)
    x_left, x_right = -0.05, 0.65
    y_low, y_high = -6.5, 0.5
    ratio = 1
    interval = 0.2
    
    ax.set_xlim([x_left, x_right])
    ax.set_ylim([y_low, y_high])
    ax.set_aspect(abs((x_right - x_left) / (y_low - y_high)) * ratio)
    
    # Labels with explicit Arial font
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    
    # Ticks and grid
    ax.set_xticks(np.arange(0 , 0.71, interval))
    ax.set_yticks(np.arange(-6, 0.1, 2))
    ax.tick_params(axis='both', labelsize=36)
    
    # Spine formatting
    for _, s in ax.spines.items():
        s.set_color('black')
        s.set_linewidth(2)
        
    # Layout adjustments
    fig.subplots_adjust(
        top=0.981,
        bottom=0.12,
        left=0.3,
        right=0.9,
        hspace=0.2,
        wspace=0.2
    )
    
    plt.savefig(out_path, dpi=300)

if __name__ == '__main__':
    file_path = 'F:/Cui_Lab/Projects/Connectional_Variability_Axis/step_03_structural_connectome_variability/sensitivity_analyses/'
    x_label = 'FC variability'
    y_label = 'Log(CMY variability)'
    
    input_file = file_path + '/yen_fc_var_sc_var_schaefer400.csv'
    out_path = file_path + '/density_plot_fc_var_sc_var_yen.png'   
    density_show(input_file,x_label,y_label,out_path)
