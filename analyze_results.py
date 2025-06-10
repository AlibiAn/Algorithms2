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

def load_tsp_coordinates(filename):
    coordinates = []
    reading_coords = False
    
    try:
        with open(filename, 'r') as file:
            for line in file:
                line = line.strip()
                
                if line == "NODE_COORD_SECTION":
                    reading_coords = True
                    continue
                
                if line == "EOF" or line.startswith("EDGE_WEIGHT_SECTION"):
                    break
                
                if reading_coords and line:
                    parts = line.split()
                    if len(parts) >= 3:
                        try:
                            x = float(parts[1])
                            y = float(parts[2])
                            coordinates.append([x, y])
                        except (ValueError, IndexError):
                            continue
        
        return np.array(coordinates)
    
    except FileNotFoundError:
        print(f"Warning: {filename} not found. Using fallback method.")
        return None
    except Exception as e:
        print(f"Error reading {filename}: {e}")
        return None

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
        hcs_data = df_agg[(df_agg['Dataset'] == dataset) & (df_agg['Algorithm'] == 'Mine')]
        hcs_data = hcs_data.sort_values('K_Clusters')
        
        axes[i].plot(hcs_data['K_Clusters'], hcs_data['Mean_Runtime_sec'], 
                    'o-', color=colors[i], linewidth=3, markersize=8, markerfacecolor='white', 
                    markeredgewidth=2, markeredgecolor=colors[i])
        axes[i].fill_between(hcs_data['K_Clusters'], 
                           hcs_data['Mean_Runtime_sec'] - hcs_data['Std_Runtime_sec'],
                           hcs_data['Mean_Runtime_sec'] + hcs_data['Std_Runtime_sec'],
                           alpha=0.2, color=colors[i])
        
        axes[i].set_title(f'{dataset}', fontsize=14, fontweight='bold')
        axes[i].set_xlabel('Number of Clusters (K)', fontsize=12)
        axes[i].set_ylabel('Runtime (seconds)', fontsize=12)
        axes[i].grid(True, alpha=0.3)
        axes[i].set_ylim(bottom=0)
    
    plt.suptitle('K-Size vs Runtime Performance for HCS Algorithm', fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'k_size_vs_runtime.png'), bbox_inches='tight')
    plt.close()

def plot_performance_vs_baseline_enhanced(df_agg):
    plt.figure(figsize=(14, 10))
    
    algorithms = ['MST-2-Approx', 'Mine']
    colors = {'MST-2-Approx': '#ff7f0e', 'Mine': '#1f77b4'}
    markers = {'MST-2-Approx': 'o', 'Mine': 's'}
    labels = {'MST-2-Approx': 'MST-2-Approx', 'Mine': 'HCS'}
    
    max_runtime = 0
    min_n = float('inf')
    max_n = 0
    
    for alg in algorithms:
        data = df_agg[df_agg['Algorithm'] == alg]
        if not data.empty:
            plt.scatter(data['Vertices'], data['Mean_Runtime_sec'], 
                       s=150, alpha=0.8, label=labels[alg], 
                       color=colors.get(alg, 'gray'), marker=markers.get(alg, 'o'),
                       edgecolors='black', linewidth=1)
            max_runtime = max(max_runtime, data['Mean_Runtime_sec'].max())
            min_n = min(min_n, data['Vertices'].min())
            max_n = max(max_n, data['Vertices'].max())
    
    n_range = np.logspace(np.log10(min_n), np.log10(max_n), 100)
    scale_factor = max_runtime / 1000
    
    n_exp = n_range[n_range <= 25]
    if len(n_exp) > 0:
        two_n = (2 ** (n_exp/8)) * scale_factor / 100
        plt.loglog(n_exp, two_n, '--', color='red', alpha=0.7, linewidth=2, label='O(2ⁿ)')
    
    n_fact_range = n_range[n_range <= 15]
    if len(n_fact_range) > 0:
        n_factorial = []
        for n in n_fact_range:
            if n <= 12:
                n_factorial.append(factorial(int(n)) * scale_factor / 1e15)
            else:
                n_factorial.append(np.nan)
        plt.loglog(n_fact_range, n_factorial, '--', color='black', alpha=0.8, linewidth=3, label='O(n!)')
    
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
    hcs_errors = []
    mst_std_errors = []
    hcs_std_errors = []
    
    for dataset in datasets:
        mst_data = df_agg[(df_agg['Dataset'] == dataset) & (df_agg['Algorithm'] == 'MST-2-Approx')]
        hcs_data = df_agg[(df_agg['Dataset'] == dataset) & (df_agg['Algorithm'] == 'Mine')]
        
        if not mst_data.empty:
            mst_errors.append(mst_data['Mean_Approx_Ratio'].iloc[0])
            mst_std_errors.append(mst_data['Std_Approx_Ratio'].iloc[0])
        else:
            mst_errors.append(0)
            mst_std_errors.append(0)
            
        if not hcs_data.empty and dataset != 'mona-lisa100K':
            best_hcs = hcs_data.loc[hcs_data['Mean_Cost'].idxmin()]
            hcs_errors.append(best_hcs['Mean_Approx_Ratio'])
            hcs_std_errors.append(best_hcs['Std_Approx_Ratio'])
        else:
            hcs_errors.append(0)
            hcs_std_errors.append(0)
    
    ax1.bar(x - width/2, mst_errors, width, yerr=mst_std_errors, capsize=5,
            label='MST-2-Approx', color='#ff7f0e', alpha=0.8)
    
    hcs_positions = []
    hcs_values = []
    hcs_errors_filtered = []
    for i, val in enumerate(hcs_errors):
        if val > 0:
            hcs_positions.append(x[i] + width/2)
            hcs_values.append(val)
            hcs_errors_filtered.append(hcs_std_errors[i])
    
    if hcs_positions:
        ax1.bar(hcs_positions, hcs_values, width, yerr=hcs_errors_filtered, capsize=5,
                label='HCS', color='#1f77b4', alpha=0.8)
    
    ax1.axhline(y=1.0, color='red', linestyle='--', alpha=0.7, linewidth=2, label='Optimal (1.0)')
    ax1.set_xlabel('Dataset', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Approximation Ratio', fontsize=12, fontweight='bold')
    ax1.set_title('Approximation Error Comparison', fontsize=14, fontweight='bold')
    ax1.set_xticks(x)
    ax1.set_xticklabels(datasets)
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    ax2.bar(x - width/2, mst_std_errors, width, label='MST-2-Approx', color='#ff7f0e', alpha=0.8)
    
    hcs_std_positions = []
    hcs_std_values = []
    for i, val in enumerate(hcs_std_errors):
        if val > 0:
            hcs_std_positions.append(x[i] + width/2)
            hcs_std_values.append(val)
    
    if hcs_std_positions:
        ax2.bar(hcs_std_positions, hcs_std_values, width, label='HCS', color='#1f77b4', alpha=0.8)
    
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
    labels = {'MST-2-Approx': 'MST-2-Approx', 'Mine': 'HCS'}
    
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
               label=labels[alg], color=colors[alg], alpha=0.8)
    
    ax.set_xlabel('Dataset', fontsize=14, fontweight='bold')
    ax.set_ylabel('Runtime (seconds)', fontsize=14, fontweight='bold')
    ax.set_title('Runtime Comparison: MST vs HCS Algorithm', fontsize=16, fontweight='bold')
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
    
    hcs_costs = []
    for dataset in datasets:
        if dataset == 'mona-lisa100K':
            hcs_costs.append(0)
        else:
            data = df_agg[(df_agg['Dataset'] == dataset) & (df_agg['Algorithm'] == 'Mine')]
            if not data.empty:
                best_data = data.loc[data['Mean_Cost'].idxmin()]
                hcs_costs.append(best_data['Mean_Cost'])
            else:
                hcs_costs.append(0)
    
    hcs_positions = []
    hcs_values = []
    for i, cost in enumerate(hcs_costs):
        if cost > 0:
            hcs_positions.append(x[i] + width)
            hcs_values.append(cost)
    
    if hcs_positions:
        ax.bar(hcs_positions, hcs_values, width, label='HCS', color='#1f77b4', alpha=0.8)
    
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
        
        if dataset != 'mona-lisa100K' and hcs_costs[i] > 0:
            hcs_pct = ((hcs_costs[i] / optimal_costs[dataset]) - 1) * 100
            ax.text(x[i] + width, hcs_costs[i], f'+{hcs_pct:.1f}%', ha='center', va='bottom', fontsize=9)
    
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'cost_comparison_with_optimal.png'), bbox_inches='tight')
    plt.close()

def plot_dataset_geometries(df_agg):
    mona_lisa_coords = load_tsp_coordinates('data/mona-lisa100K.tsp')
    kz9976_coords = load_tsp_coordinates('data/kz9976.tsp')
    
    if mona_lisa_coords is None:
        mona_lisa_coords = load_tsp_coordinates('data/monalisa.tsp')
    if mona_lisa_coords is None:
        mona_lisa_coords = load_tsp_coordinates('data/mona_lisa.tsp')
    
    if kz9976_coords is None:
        kz9976_coords = load_tsp_coordinates('data/kz.tsp')
    
    fig, axes = plt.subplots(1, 2, figsize=(16, 8))
    
    if mona_lisa_coords is not None:
        axes[0].scatter(mona_lisa_coords[:, 0], mona_lisa_coords[:, 1], 
                       s=0.3, alpha=0.4, color='#1f77b4', rasterized=True)
        axes[0].set_title('mona-lisa100K Cities (Actual Data)', fontsize=14, fontweight='bold')
        print(f"Loaded {len(mona_lisa_coords)} mona-lisa coordinates")
    else:
        np.random.seed(42)
        clusters = 12
        points_per_cluster = 8333
        mona_lisa_coords = []
        
        for i in range(clusters):
            if i < 4:
                center = [400 + (i-1.5)*100, 500 + (i%2)*50]
                spread = 30
            elif i < 8:
                center = [300 + (i-4)*100, 600 + ((i-4)%2)*40]
                spread = 45
            else:
                center = [200 + (i-8)*200, 400 + ((i-8)%2)*100]
                spread = 60
            
            cluster_points = np.array(center) + np.random.randn(points_per_cluster, 2) * spread
            mona_lisa_coords.append(cluster_points)
        
        mona_lisa_coords = np.vstack(mona_lisa_coords)
        axes[0].scatter(mona_lisa_coords[:, 0], mona_lisa_coords[:, 1], 
                       s=0.3, alpha=0.4, color='#1f77b4', rasterized=True)
        axes[0].set_title('mona-lisa100K Cities (Simulated)', fontsize=14, fontweight='bold')
        print("Using simulated mona-lisa coordinates (TSP file not found)")
    
    axes[0].set_xlabel('X Coordinate', fontsize=12)
    axes[0].set_ylabel('Y Coordinate', fontsize=12)
    axes[0].grid(True, alpha=0.3)
    
    if kz9976_coords is not None:
        axes[1].scatter(kz9976_coords[:, 0], kz9976_coords[:, 1], 
                       s=0.5, alpha=0.5, color='#ff7f0e', rasterized=True)
        axes[1].set_title('kz9976 Cities (Actual Data)', fontsize=14, fontweight='bold')
        print(f"Loaded {len(kz9976_coords)} kz9976 coordinates")
    else:
        np.random.seed(123)
        t = np.linspace(0, 12*np.pi, 9976)
        x = t * np.cos(t/4) + np.random.randn(9976) * 200
        y = t * np.sin(t/4) + np.random.randn(9976) * 150
        for i in range(0, len(x), 1000):
            end_idx = min(i+1000, len(x))
            cluster_center = [np.random.randn()*500, np.random.randn()*500]
            x[i:end_idx] += cluster_center[0]
            y[i:end_idx] += cluster_center[1]
        
        kz9976_coords = np.vstack((x, y)).T
        axes[1].scatter(kz9976_coords[:, 0], kz9976_coords[:, 1], 
                       s=0.5, alpha=0.5, color='#ff7f0e', rasterized=True)
        axes[1].set_title('kz9976 Cities (Simulated)', fontsize=14, fontweight='bold')
        print("Using simulated kz9976 coordinates (TSP file not found)")
    
    axes[1].set_xlabel('X Coordinate', fontsize=12)
    axes[1].set_ylabel('Y Coordinate', fontsize=12)
    axes[1].grid(True, alpha=0.3)
    
    fig.suptitle("Datasets",
                fontsize=11, y=0.02, wrap=True)
    
    plt.tight_layout(rect=[0, 0.08, 1, 0.96])
    plt.savefig(os.path.join(OUTPUT_DIR, 'dataset_geometries_comparison.png'), 
               bbox_inches='tight', dpi=300)
    plt.close()
    
    print("6. dataset_geometries_comparison.png - Scatter plots comparing dataset geometries")

def main():
    setup_environment()
    df, df_agg = load_data()
    
    print("Generating all 6 runtime analysis plots...")
    
    plot_k_size_vs_runtime(df_agg)
    plot_performance_vs_baseline_enhanced(df_agg)
    plot_error_std_comparison_fixed(df_agg)
    plot_runtime_comparison_all(df_agg)
    plot_cost_comparison_with_optimal(df_agg)
    plot_dataset_geometries(df_agg)
    
    print(f"All 6 plots saved in '{OUTPUT_DIR}' directory!")
    print("Complete set of analysis plots:")
    print("1. k_size_vs_runtime.png - K-size vs runtime with connected lines")
    print("2. performance_vs_baseline.png - Enhanced with n! and n! log n baselines")
    print("3. error_std_comparison.png - FIXED: MST-2-Approx bars now showing")
    print("4. runtime_comparison_all.png - Runtime comparison in seconds")
    print("5. cost_comparison_with_optimal.png - Cost comparison with optimal solutions")
    print("6. dataset_geometries_comparison.png - Dataset geometry comparison")

if __name__ == '__main__':
    main()
