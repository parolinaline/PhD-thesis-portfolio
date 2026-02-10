import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

# Create the data from your heatmap
data = {
    'comparison': [
        'CP4 vs CP10', 'CP4 vs CP10', 'CP4 vs CP10', 'CP4 vs CP10', 'CP4 vs CP10', 'CP4 vs CP10', 'CP4 vs CP10', 'CP4 vs CP10', 'CP4 vs CP10', 'CP4 vs CP10', 'CP4 vs CP10',
        'CP5 vs CP10', 'CP5 vs CP10', 'CP5 vs CP10', 'CP5 vs CP10', 'CP5 vs CP10', 'CP5 vs CP10', 'CP5 vs CP10', 'CP5 vs CP10', 'CP5 vs CP10', 'CP5 vs CP10', 'CP5 vs CP10',
        'CP5 vs CP4', 'CP5 vs CP4', 'CP5 vs CP4', 'CP5 vs CP4', 'CP5 vs CP4', 'CP5 vs CP4', 'CP5 vs CP4', 'CP5 vs CP4', 'CP5 vs CP4', 'CP5 vs CP4', 'CP5 vs CP4',
        'CP4 vs C.j. 81-176', 'CP4 vs C.j. 81-176', 'CP4 vs C.j. 81-176', 'CP4 vs C.j. 81-176', 'CP4 vs C.j. 81-176', 'CP4 vs C.j. 81-176', 'CP4 vs C.j. 81-176', 'CP4 vs C.j. 81-176', 'CP4 vs C.j. 81-176', 'CP4 vs C.j. 81-176', 'CP4 vs C.j. 81-176',
        'CP5 vs C.j. 81-176', 'CP5 vs C.j. 81-176', 'CP5 vs C.j. 81-176', 'CP5 vs C.j. 81-176', 'CP5 vs C.j. 81-176', 'CP5 vs C.j. 81-176', 'CP5 vs C.j. 81-176', 'CP5 vs C.j. 81-176', 'CP5 vs C.j. 81-176', 'CP5 vs C.j. 81-176', 'CP5 vs C.j. 81-176',
        'CP10 vs C.j. 81-176', 'CP10 vs C.j. 81-176', 'CP10 vs C.j. 81-176', 'CP10 vs C.j. 81-176', 'CP10 vs C.j. 81-176', 'CP10 vs C.j. 81-176', 'CP10 vs C.j. 81-176', 'CP10 vs C.j. 81-176', 'CP10 vs C.j. 81-176', 'CP10 vs C.j. 81-176', 'CP10 vs C.j. 81-176',
        'CP4 vs Agro', 'CP4 vs Agro', 'CP4 vs Agro', 'CP4 vs Agro', 'CP4 vs Agro', 'CP4 vs Agro', 'CP4 vs Agro', 'CP4 vs Agro', 'CP4 vs Agro', 'CP4 vs Agro', 'CP4 vs Agro',
        'CP5 vs Agro', 'CP5 vs Agro', 'CP5 vs Agro', 'CP5 vs Agro', 'CP5 vs Agro', 'CP5 vs Agro', 'CP5 vs Agro', 'CP5 vs Agro', 'CP5 vs Agro', 'CP5 vs Agro', 'CP5 vs Agro',
        'CP10 vs Agro', 'CP10 vs Agro', 'CP10 vs Agro', 'CP10 vs Agro', 'CP10 vs Agro', 'CP10 vs Agro', 'CP10 vs Agro', 'CP10 vs Agro', 'CP10 vs Agro', 'CP10 vs Agro', 'CP10 vs Agro'
    ],
    'gene': ['virB2', 'virB3', 'virB4', 'virB5', 'virB6', 'virB7', 'virB8', 'virB9', 'virB10', 'virB11', 'virD'] * 9,
    'identity': [
        # CP4 vs CP10
        36.53, np.nan, 23.45, 23.52, 23.75, np.nan, 24.40, 28.91, 29.41, 29.52, 24.29,
        # CP5 vs CP10
        0, 0, 23.49, 27.20, np.nan, 30, 0, 40, 31.07, 30.54, 24.50,
        # CP5 vs CP4
        42.85, np.nan, 23.24, 56.25, np.nan, np.nan, 22.56, 21.47, 27.52, 26.07, 25.13,
        # CP4 vs C. jejuni 81-176
        82.75, np.nan, 97.93, 93.86, 93.30, np.nan, 100, 100, 100, 100, 100,
        # CP5 vs C. jejuni 81-176
        0, 28.12, 23.35, 28.07, np.nan, 30, 23.49, 21.47, 27.52, 26.07, 25.13,
        # CP10 vs C. jejuni 81-176
        35.84, 45, 22.75, 91.13, 24.68, 100, 24.43, 28.91, 29.41, 29.52, 24.29,
        # CP4 vs Agrobacterium
        0, np.nan, 21.72, 0, 0, np.nan, 0.00, 24.08, 30.10, 25.78, 27.39,
        # CP5 vs Agrobacterium
        46.25, 44.68, 44.98, 24.51, np.nan, 24.44, 31.60, 31.60, 51.11, 49.35, 34.75,
        # CP10 vs Agrobacterium
        37.93, 0, 21.07, 40, 25.49, 0, 25.35, 23.07, 25.78, 30.14, 23.34
    ]
}

df = pd.DataFrame(data)

# Create a column to distinguish between missing genes (NaN) and 0% identity
df['data_type'] = df['identity'].apply(lambda x: 'Missing gene' if pd.isna(x) else
                                        ('No similarity' if x == 0 else 'Present'))

# Add category for coloring
def get_category(comp):
    if 'CP4 vs CP' in comp or 'CP5 vs CP' in comp:
        return 'Local plasmids'
    elif 'vs C.j.' in comp:
        return 'vs C.j. 81-176'
    else:
        return 'vs Agrobacterium'

df['category'] = df['comparison'].apply(get_category)

# Create faceted plot
fig, axes = plt.subplots(3, 3, figsize=(15, 12))
axes = axes.flatten()

# Reorder comparisons: local plasmids first, then reference strain, then Agrobacterium
comparisons_ordered = [
    'CP5 vs CP4', 'CP5 vs CP10', 'CP4 vs CP10',  # Local plasmids
    'CP4 vs C.j. 81-176', 'CP5 vs C.j. 81-176', 'CP10 vs C.j. 81-176',  # Reference strain
    'CP4 vs Agro', 'CP5 vs Agro', 'CP10 vs Agro'  # Agrobacterium
]
category_colors = {'Local plasmids': '#1b9e77', 'vs C.j. 81-176': '#d95f02', 'vs Agrobacterium': '#7570b3'}

for idx, comparison in enumerate(comparisons_ordered):
    ax = axes[idx]
    data_subset = df[df['comparison'] == comparison].copy()

    # Separate data with values from missing data
    data_with_values = data_subset[data_subset['data_type'] != 'Missing gene'].copy()
    data_missing = data_subset[data_subset['data_type'] == 'Missing gene'].copy()

    # Sort by identity for better visualization
    data_with_values = data_with_values.sort_values('identity', ascending=False)

    # Plot data with identity values
    color = data_with_values['category'].iloc[0] if len(data_with_values) > 0 else 'gray'

    if len(data_with_values) > 0:
        ax.scatter(data_with_values['identity'], range(len(data_with_values)),
                  s=120, alpha=0.8, color=category_colors[color],
                  edgecolors='black', linewidth=1, label='Gene present', zorder=3)

    # Plot missing genes as X marks
    if len(data_missing) > 0:
        ax.scatter([102] * len(data_missing), range(len(data_with_values),
                   len(data_with_values) + len(data_missing)),
                  s=200, marker='x', color='red', linewidth=2.5,
                  label='Gene absent', zorder=3)

    # Combine gene labels
    all_genes = list(data_with_values['gene']) + list(data_missing['gene'])
    ax.set_yticks(range(len(all_genes)))
    ax.set_yticklabels(all_genes, fontsize=9)

    # Labels and formatting
    ax.set_xlabel('Identity (%)', fontsize=9)
    ax.set_title(comparison, fontsize=10, fontweight='bold')
    ax.set_xlim(-5, 110)
    ax.set_ylim(-0.5, len(all_genes) - 0.5)
    ax.grid(axis='x', alpha=0.3)

    # Add vertical line at 50% for reference
    ax.axvline(50, color='gray', linestyle='--', alpha=0.5, linewidth=1)

    # Add legend only to first subplot
    if idx == 0:
        ax.legend(loc='lower right', fontsize=8)

# Remove extra subplots
for idx in range(len(comparisons_ordered), len(axes)):
    fig.delaxes(axes[idx])

# Create main legend
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
legend_elements = [
    Patch(facecolor=category_colors['Local plasmids'], edgecolor='black', label='Local plasmids'),
    Patch(facecolor=category_colors['vs C.j. 81-176'], edgecolor='black', label='vs C.j. 81-176'),
    Patch(facecolor=category_colors['vs Agrobacterium'], edgecolor='black', label='vs Agrobacterium'),
    Line2D([0], [0], marker='x', color='w', markerfacecolor='red', markersize=10,
           markeredgewidth=2.5, label='Gene absent (NA)')
]

fig.legend(handles=legend_elements, loc='upper center', bbox_to_anchor=(0.5, 0.98),
          ncol=4, frameon=True, fontsize=10)

plt.suptitle('T4SS Gene Identity Across Plasmid Comparisons', fontsize=14, fontweight='bold', y=0.995)
plt.tight_layout(rect=[0, 0, 1, 0.96])

# Save
plt.savefig('t4ss_faceted_dotplot.png', dpi=300, bbox_inches='tight')
plt.savefig('t4ss_faceted_dotplot.pdf', dpi=300, bbox_inches='tight')

print("Faceted dot plot created successfully!")
print(f"Total data points: {len(df)}")
print(f"Missing genes (NA): {df['data_type'].value_counts().get('Missing gene', 0)}")
print(f"Genes with no similarity (0%): {df['data_type'].value_counts().get('No similarity', 0)}")
print(f"Genes with identity data: {df['data_type'].value_counts().get('Present', 0)}")
