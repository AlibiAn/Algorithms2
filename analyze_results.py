import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import os
from math import factorial

CSV_FILE = 'tsp_results.csv'
OUTPUT_DIR = 'runtime_analysis_plots'

def setup_environment():
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
    sns.set_style('whitegrid')
    plt.rcParams['figure.facecolor'] = 'white'
    plt.rcParams['savefig.dpi'] = 300
    plt.rcParams['font.size'] = 12

def load_data():
    df = pd.read_csv(CSV_FILE)
    df_agg = df.groupby(['Dataset', 'Algorithm', 'Vertices', 'K_Clusters']).agg({
        'Tour_Cost': ['mean', 'std'],
        'Runtime_ms': ['mean', 'std'],
        'Approximation_Ratio': ['mean', 'std']
    }).reset_index()
    df_agg.columns = ['Dataset', 'Algorithm', 'Vertices', 'K_Clusters', 'Mean_Cost', 'Std_Cost', 'Mean_Runtime_ms', 'Std_Runtime_ms', 'Mean_Approx_Ratio', 'Std_Approx_Ratio']
    df_agg['Std_Cost'] = df_agg['Std_Cost'].fillna(0)
    df_agg['Std_Runtime_ms'] = df_agg['Std_Runtime_ms'].fillna(0)
    df_agg['Mean_Runtime_sec'] = df_agg['Mean_Runtime_ms'] / 1000
    df_agg['Std_Runtime_sec'] = df_agg['Std_Runtime_ms'] / 1000
    return df, df_agg

def plot_k_size_vs_runtime(df_agg):
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    datasets = ['a280', 'xql662', 'kz9976']
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c']
    
    for i, dataset in enumerate(datasets):
        mine_data = df_agg[(df_agg['Dataset'] == dataset) & (df_agg['Algorithm'] == 'Mine')]
        mine_data = mine_data.sort_values('K_Clusters')
        
        axes[i].plot(mine_data['K_Clusters'], mine_data['Mean_Runtime_sec'], 
                    'o-', color=colors[i], linewidth=3, markersize=8, markerfacecolor='white', 
                    markeredgewidth=2, markeredgecolor=colors[i])
        axes[i].fill_between(mine_data['K_Clusters'], 
                           mine_data['Mean_Runtime_sec'] - mine_data['Std_Runtime_sec'],
                           mine_data['Mean_Runtime_sec'] + mine_data['Std_Runtime_sec'],
                           alpha=0.2, color=colors[i])
        
        axes[i].set_title(f'{dataset}', fontsize=14, fontweight='bold')
        axes[i].set_xlabel('Number of Clusters (K)', fontsize=12)
        axes[i].set_ylabel('Runtime (seconds)', fontsize=12)
        axes[i].grid(True, alpha=0.3)
        axes[i].set_ylim(bottom=0)
    
    plt.suptitle('K-Size vs Runtime Performance for My Algorithm', fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'k_size_vs_runtime.png'), bbox_inches='tight')
    plt.close()

def plot_performance_vs_baseline_enhanced(df_agg):
    plt.figure(figsize=(14, 10))
    
    algorithms = ['MST-2-Approx', 'Mine']
    colors = {'MST-2-Approx': '#ff7f0e', 'Mine': '#1f77b4'}
    markers = {'MST-2-Approx': 'o', 'Mine': 's'}
    
    max_runtime = 0
    min_n = float('inf')
    max_n = 0
    
    for alg in algorithms:
        data = df_agg[df_agg['Algorithm'] == alg]
        if not data.empty:
            plt.scatter(data['Vertices'], data['Mean_Runtime_sec'], 
                       s=150, alpha=0.8, label=f'{alg}', 
                       color=colors.get(alg, 'gray'), marker=markers.get(alg, 'o'),
                       edgecolors='black', linewidth=1)
            max_runtime = max(max_runtime, data['Mean_Runtime_sec'].max())
            min_n = min(min_n, data['Vertices'].min())
            max_n = max(max_n, data['Vertices'].max())
    
    n_range = np.logspace(np.log10(min_n), np.log10(max_n), 100)
    scale_factor = max_runtime / 1000
    
    # Exponential baseline
    n_exp = n_range[n_range <= 25]
    if len(n_exp) > 0:
        two_n = (2 ** (n_exp/8)) * scale_factor / 100
        plt.loglog(n_exp, two_n, '--', color='red', alpha=0.7, linewidth=2, label='O(2ⁿ)')
    
    # NON-POLYNOMIAL BASELINES - FACTORIAL
    n_fact_range = n_range[n_range <= 15]
    if len(n_fact_range) > 0:
        n_factorial = []
        for n in n_fact_range:
            if n <= 12:
                n_factorial.append(factorial(int(n)) * scale_factor / 1e15)
            else:
                n_factorial.append(np.nan)
        plt.loglog(n_fact_range, n_factorial, '--', color='black', alpha=0.8, linewidth=3, label='O(n!)')
    
    # NON-POLYNOMIAL BASELINES - FACTORIAL LOG N
    if len(n_fact_range) > 0:
        n_fact_log_n = []
        for n in n_fact_range:
            if n <= 10:
                n_fact_log_n.append(factorial(int(n)) * np.log2(n) * scale_factor / 1e16)
            else:
                n_fact_log_n.append(np.nan)
        plt.loglog(n_fact_range, n_fact_log_n, '--', color='brown', alpha=0.8, linewidth=3, label='O(n! log n)')
    
    plt.xlabel('Dataset Size (vertices)', fontsize=14, fontweight='bold')
    plt.ylabel('Runtime (seconds)', fontsize=14, fontweight='bold')
    plt.title('Algorithm Performance vs Extended Theoretical Complexity Baselines', fontsize=16, fontweight='bold')
    plt.legend(fontsize=10, loc='upper left')
    plt.grid(True, which='both', alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'performance_vs_baseline.png'), bbox_inches='tight')
    plt.close()

def plot_error_std_comparison_fixed(df_agg):
    datasets = ['a280', 'xql662', 'kz9976', 'mona-lisa100K']
    algorithms = ['MST-2-Approx', 'Mine']
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    
    x = np.arange(len(datasets))
    width = 0.35
    
    mst_errors = []
    mine_errors = []
    mst_std_errors = []
    mine_std_errors = []
    
    for dataset in datasets:
        mst_data = df_agg[(df_agg['Dataset'] == dataset) & (df_agg['Algorithm'] == 'MST-2-Approx')]
        mine_data = df_agg[(df_agg['Dataset'] == dataset) & (df_agg['Algorithm'] == 'Mine')]
        
        # CRITICAL BUG FIX: Added [0] to .iloc
        if not mst_data.empty:
            mst_errors.append(mst_data['Mean_Approx_Ratio'].iloc[0])
            mst_std_errors.append(mst_data['Std_Approx_Ratio'].iloc[0])
        else:
            mst_errors.append(0)
            mst_std_errors.append(0)
            
        if not mine_data.empty and dataset != 'mona-lisa100K':
            best_mine = mine_data.loc[mine_data['Mean_Cost'].idxmin()]
            mine_errors.append(best_mine['Mean_Approx_Ratio'])
            mine_std_errors.append(best_mine['Std_Approx_Ratio'])
        else:
            mine_errors.append(0)
            mine_std_errors.append(0)
    
    # Plot approximation ratios
    ax1.bar(x - width/2, mst_errors, width, yerr=mst_std_errors, capsize=5,
            label='MST-2-Approx', color='#ff7f0e', alpha=0.8)
    
    mine_positions = []
    mine_values = []
    mine_errors_filtered = []
    for i, val in enumerate(mine_errors):
        if val > 0:
            mine_positions.append(x[i] + width/2)
            mine_values.append(val)
            mine_errors_filtered.append(mine_std_errors[i])
    
    if mine_positions:
        ax1.bar(mine_positions, mine_values, width, yerr=mine_errors_filtered, capsize=5,
                label='My Algorithm', color='#1f77b4', alpha=0.8)
    
    ax1.axhline(y=1.0, color='red', linestyle='--', alpha=0.7, linewidth=2, label='Optimal (1.0)')
    ax1.set_xlabel('Dataset', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Approximation Ratio', fontsize=12, fontweight='bold')
    ax1.set_title('Approximation Error Comparison', fontsize=14, fontweight='bold')
    ax1.set_xticks(x)
    ax1.set_xticklabels(datasets)
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # FIXED: Standard deviation comparison - MST bars should now show
    ax2.bar(x - width/2, mst_std_errors, width, label='MST-2-Approx', color='#ff7f0e', alpha=0.8)
    
    mine_std_positions = []
    mine_std_values = []
    for i, val in enumerate(mine_std_errors):
        if val > 0:
            mine_std_positions.append(x[i] + width/2)
            mine_std_values.append(val)
    
    if mine_std_positions:
        ax2.bar(mine_std_positions, mine_std_values, width, label='My Algorithm', color='#1f77b4', alpha=0.8)
    
    ax2.set_xlabel('Dataset', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Standard Deviation', fontsize=12, fontweight='bold')
    ax2.set_title('Solution Consistency (Lower = More Consistent)', fontsize=14, fontweight='bold')
    ax2.set_xticks(x)
    ax2.set_xticklabels(datasets)
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'error_std_comparison.png'), bbox_inches='tight')
    plt.close()

def plot_runtime_comparison_all(df_agg):
    datasets = ['a280', 'xql662', 'kz9976', 'mona-lisa100K']
    algorithms = ['MST-2-Approx', 'Mine']
    colors = {'MST-2-Approx': '#ff7f0e', 'Mine': '#1f77b4'}
    
    fig, ax = plt.subplots(figsize=(16, 10))
    
    x = np.arange(len(datasets))
    width = 0.35
    
    for i, alg in enumerate(algorithms):
        runtimes = []
        std_runtimes = []
        
        for dataset in datasets:
            data = df_agg[(df_agg['Dataset'] == dataset) & (df_agg['Algorithm'] == alg)]
            
            if not data.empty:
                if alg == 'Mine':
                    best_data = data.loc[data['Mean_Cost'].idxmin()]
                    runtimes.append(best_data['Mean_Runtime_sec'])
                    std_runtimes.append(best_data['Std_Runtime_sec'])
                else:
                    runtimes.append(data['Mean_Runtime_sec'].iloc[0])
                    std_runtimes.append(data['Std_Runtime_sec'].iloc[0])
            else:
                runtimes.append(0)
                std_runtimes.append(0)
        
        if alg == 'Mine':
            runtimes[-1] = 0
            std_runtimes[-1] = 0
        
        mask = np.array(runtimes) > 0
        positions = x[mask] + i * width - width/2
        values = np.array(runtimes)[mask]
        errors = np.array(std_runtimes)[mask]
        
        ax.bar(positions, values, width, yerr=errors, capsize=5,
               label=alg, color=colors[alg], alpha=0.8)
    
    ax.set_xlabel('Dataset', fontsize=14, fontweight='bold')
    ax.set_ylabel('Runtime (seconds)', fontsize=14, fontweight='bold')
    ax.set_title('Runtime Comparison: MST vs My Algorithm', fontsize=16, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(datasets)
    ax.set_yscale('log')
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3, which='both')
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'runtime_comparison_all.png'), bbox_inches='tight')
    plt.close()

def plot_cost_comparison_with_optimal(df_agg):
    datasets = ['a280', 'xql662', 'kz9976', 'mona-lisa100K']
    optimal_costs = {
        'a280': 2579.0,
        'xql662': 2513.0,
        'kz9976': 1061882.0,
        'mona-lisa100K': 5757191.0
    }
    
    fig, ax = plt.subplots(figsize=(16, 10))
    
    x = np.arange(len(datasets))
    width = 0.25
    
    opt_costs = [optimal_costs[d] for d in datasets]
    ax.bar(x - width, opt_costs, width, label='Optimal', color='green', alpha=0.8)
    
    mst_costs = []
    for dataset in datasets:
        data = df_agg[(df_agg['Dataset'] == dataset) & (df_agg['Algorithm'] == 'MST-2-Approx')]
        if not data.empty:
            mst_costs.append(data['Mean_Cost'].iloc[0])
        else:
            mst_costs.append(0)
    
    ax.bar(x, mst_costs, width, label='MST-2-Approx', color='#ff7f0e', alpha=0.8)
    
    mine_costs = []
    for dataset in datasets:
        if dataset == 'mona-lisa100K':
            mine_costs.append(0)
        else:
            data = df_agg[(df_agg['Dataset'] == dataset) & (df_agg['Algorithm'] == 'Mine')]
            if not data.empty:
                best_data = data.loc[data['Mean_Cost'].idxmin()]
                mine_costs.append(best_data['Mean_Cost'])
            else:
                mine_costs.append(0)
    
    mine_positions = []
    mine_values = []
    for i, cost in enumerate(mine_costs):
        if cost > 0:
            mine_positions.append(x[i] + width)
            mine_values.append(cost)
    
    if mine_positions:
        ax.bar(mine_positions, mine_values, width, label='My Algorithm', color='#1f77b4', alpha=0.8)
    
    ax.set_xlabel('Dataset', fontsize=14, fontweight='bold')
    ax.set_ylabel('Tour Cost', fontsize=14, fontweight='bold')
    ax.set_title('Tour Cost Comparison with Optimal Solutions', fontsize=16, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(datasets)
    ax.set_yscale('log')
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3, which='both')
    
    for i, dataset in enumerate(datasets):
        if mst_costs[i] > 0:
            mst_pct = ((mst_costs[i] / optimal_costs[dataset]) - 1) * 100
            ax.text(x[i], mst_costs[i], f'+{mst_pct:.1f}%', ha='center', va='bottom', fontsize=9)
        
        if dataset != 'mona-lisa100K' and mine_costs[i] > 0:
            mine_pct = ((mine_costs[i] / optimal_costs[dataset]) - 1) * 100
            ax.text(x[i] + width, mine_costs[i], f'+{mine_pct:.1f}%', ha='center', va='bottom', fontsize=9)
    
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'cost_comparison_with_optimal.png'), bbox_inches='tight')
    plt.close()

def main():
    setup_environment()
    df, df_agg = load_data()
    
    print("Generating all 5 runtime analysis plots...")
    
    plot_k_size_vs_runtime(df_agg)
    plot_performance_vs_baseline_enhanced(df_agg)
    plot_error_std_comparison_fixed(df_agg)
    plot_runtime_comparison_all(df_agg)
    plot_cost_comparison_with_optimal(df_agg)
    
    print(f"All 5 plots saved in '{OUTPUT_DIR}' directory!")
    print("Complete set of analysis plots:")
    print("1. k_size_vs_runtime.png - K-size vs runtime with connected lines")
    print("2. performance_vs_baseline.png - Enhanced with n! and n! log n baselines")
    print("3. error_std_comparison.png - FIXED: MST-2-Approx bars now showing")
    print("4. runtime_comparison_all.png - Runtime comparison in seconds")
    print("5. cost_comparison_with_optimal.png - Cost comparison with optimal solutions")

if __name__ == '__main__':
    main()
