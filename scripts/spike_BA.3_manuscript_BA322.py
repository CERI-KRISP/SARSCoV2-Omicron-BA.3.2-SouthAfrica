import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.lines import Line2D
from matplotlib.colors import Normalize

#reviewed

# mutation_positions = [
#     2, 28, 31, 67, 69, 95, 101, 136, 157, 164, 172, 187,
#     211, 212, 214, 242, 251, 326, 339, 348, 356, 371, 373, 375,
#     403, 405, 408, 417, 435, 440, 445, 446, 452, 460, 477,
#     478, 484, 496, 498, 501, 529, 554, 575, 583, 614, 625,
#     641, 642, 654, 655, 679, 681, 688, 704, 764, 795, 796,
#     852, 939, 954, 969, 1184
# ]

# mutation_names = [
#     'P9L', 'R21T', 'P26L', 'A67V', 'del69-70', 'T95I', 'I101T',
#     'del136-147', 'F157S', 'N164K', 'S172F', 'K187T',
#     'del211', 'L212I', 'Ins214-ASDT', 'del242-243', 'P251S', 'I326V', 'G339Y',
#     'A348P', 'K356T', 'S371F', 'S373P', 'S375F', 'R403K',
#     'D405N', 'R408S', 'K417N', 'A435S', 'N440R', 'V445A',
#     'G446D', 'L452W', 'N460K', 'S477N', 'T478N', 'E484K',
#     'G496S', 'Q498R', 'N501Y', 'K529N', 'E554D', 'A575S',
#     'E583D', 'D614G', 'H625R', 'N641K', 'V642G', 'E654K',
#     'H655Y', 'N679R', 'P681H', 'A688D', 'S704L', 'N764K',
#     'K795T', 'D796Y', 'A852K', 'S939F', 'Q954H', 'N969K',
#     'D1184E'
# ]

###UPDATED

mutation_positions = [
    9, 21, 26, 67, 69, 95, 101, 136, 157, 164, 172, 187,
    211, 212, 214, 242, 251, 326, 339, 348, 356, 371, 373, 375,
    403, 405, 408, 417, 435, 440, 445, 446, 452, 460, 477,
    478, 484, 493, 496, 498, 501, 505, 529, 554, 575, 583, 614, 625,
    641, 642, 654, 655, 679, 681, 688, 704, 764, 795, 796,
    852, 939, 954, 969, 1184
]

mutation_names = [
    'P9L', 'R21T', 'P26L', 'A67V', 'del69-70', 'T95I', 'I101T',
    'del136-147', 'F157S', 'N164K', 'S172F', 'K187T',
    'del211', 'L212I', 'Ins214-ASDT', 'del242-243', 'P251S', 'I326V', 'G339Y/D339Y',
    'A348P', '*K356T', 'S371F', 'S373P', 'S375F', 'R403K',
    'D405N', 'R408S', 'K417N', 'A435S', 'N440R/K440R', 'V445A',
    'G446D/S446D', 'L452W', 'N460K', 'S477N', 'T478N/K478N', 'E484K', 'R493Q',
    'G496S', 'Q498R', 'N501Y', 'H505Y', 'K529N', 'E554D', '*A575S',
    'E583D', 'D614G', 'H625R', 'N641K', 'V642G', 'E654K',
    'H655Y', 'N679R/K679R', '*P681H', 'A688D', 'S704L', 'N764K',
    'K795T', 'D796Y', 'A852K', 'S939F', 'Q954H', 'N969K',
    'D1184E'
]


# Updated regions
regions = {
    'SP': (1, 13), 
    'NTD': (14, 306),
    '': (307, 333),
    'RBD': (334, 528),
    'SD1/SD2': (529, 685),
    'S2': (686, 1211),
    '.': (1212, 1273)
}

# Full names of regions
region_full_names = {
    'SP': 'signal peptide',
    'NTD': 'N-Terminal Domain (NTD)',
    'RBD': 'Receptor-Binding Domain (RBD)',
    'SD1/SD2': 'Subdomain 1/Subdomain 2 (SD1/SD2)',
    'S2': 'Subunit 2 (S2)'
}

# Color for each region
# region_colors = {
#     'SP': 'red',
#     'NTD': 'cyan',
#     'RBD': 'lightgreen',
#     'SD1/SD2': 'orange',
#     'S2': 'yellow'
# }

region_colors = {
    'SP': '#d94d4d',      # red tone
    'NTD': '#4da6d9',     # blue tone
    'RBD': '#4dd97b',     # green tone
    'SD1/SD2': '#d94dc7', # magenta tone
    'S2': '#d9cf4d'       # yellow tone
}

# Default color for regions not in the list
default_color = 'grey'  # Change to the desired default color

# Plotting
fig, ax = plt.subplots(figsize=(15, 8))  # Adjust height as needed

# Define a normalization instance for coloring the regions
norm = Normalize(vmin=0, vmax=1)

# Iterate over the entire gene range and apply default color to unlabeled areas
gene_start = 0
gene_end = 1274
for region, (start, end) in regions.items():
    ax.add_patch(patches.Rectangle((gene_start, 0), start - gene_start, 1, edgecolor='black', linewidth=0,
                                   facecolor=default_color, alpha=0.2, transform=ax.transData))
    gene_start = end

# Add grey color to the end of the gene range not covered by labeled regions
last_region_end = max([end for _, end in regions.values()])
ax.add_patch(patches.Rectangle((last_region_end, 0), gene_end - last_region_end, 1, edgecolor=None, linewidth=0,
                               facecolor=default_color, alpha=0.2, transform=ax.transData))

# Add black border around each labeled region
for region, (start, end) in regions.items():
    height = 1
    ax.add_patch(patches.Rectangle((start, 0), end - start, height, edgecolor='black', linewidth=1,
                                   facecolor=region_colors.get(region, default_color), alpha=0.5, transform=ax.transData))

# Plotting region labels
for region, (start, end) in regions.items():
    ax.text((start + end) / 2, 0.5, region, ha='center', va='center', color='black')

# Plotting mutation positions
ax.scatter(mutation_positions, [0.008]*len(mutation_positions), marker='o', color='red', s=8)

highlight_mutations = {   'P9L', 'R21T', 'P26L', 'I101T',
    'del136-147', 'F157S', 'N164K', 'S172F', 'K187T',
    'Ins214-ASDT', 'del242-243', 'P251S', 'I326V', 'G339Y/D339Y',
    'A348P', '*K356T', 'R403K',
    'R408S', 'A435S', 'N440R/K440R', 'V445A',
    'G446D/S446D', 'L452W', 'N460K', 'T478N/K478N', 'R493Q',
    'G496S','H505Y', 'K529N', 'E554D', '*A575S',
    'E583D', 'H625R', 'N641K', 'V642G', 'E654K',
    'N679R/K679R', 'A688D', 'S704L',
    'K795T', 'A852K', 'S939F',
    'D1184E'}


# Adding mutation labels with connecting lines
label_y_dict = {}  # Dictionary to store label positions by region
for x, label in zip(mutation_positions, mutation_names):
    for region, (start, end) in regions.items():
        if start <= x <= end:
            if region not in label_y_dict:
                label_y_dict[region] = -1.8  # Initialize label position for new region
            label_y = label_y_dict[region] - 0.05  # Adjust label position vertically
            ax.text(x, label_y, label, ha='right', va='top', color='red' if label in highlight_mutations else 'black', fontsize=8)  # Set va to 'top'
            ax.plot([x, x], [0, label_y], color='black', linestyle='-', linewidth=0.4, transform=ax.transData)  # Connect dot to label
            label_y_dict[region] -= 0.3  # Decrement label y-coordinate by 0.1 for the next label

# Add horizontal lines and annotations for S1 and S2
ax.axhline(y=3.5, xmin=13/1273, xmax=685/1273, color='lightgreen', linestyle='-', linewidth=3)
ax.text(350, 3.6, 'S1', ha='center', va='bottom', color='black', fontsize=10)
ax.axhline(y=3.5, xmin=686/1273, xmax=1211/1273, color='green', linestyle='-', linewidth=3)
ax.text(1000, 3.6, 'S2', ha='center', va='bottom', color='black', fontsize=10)
ax.axhline(y=4.5, xmin=0/1273, xmax=1273/1273, color='green', linestyle='-', linewidth=3)
ax.text(636, 4.6, 'BA.3.2.2 spike (aa)', ha='center', va='bottom', color='black', fontsize=10)


# Add a dotted line and label for cleavage site
ax.plot([685, 685], [0, 2], color='black', linestyle=':', linewidth=1)
ax.text(685, 2.1, 'S1/S2 cleavage site 685-686', ha='center', va='bottom', color='black', fontsize=10)
ax.plot([815, 815], [0, 2.6], color='black', linestyle=':', linewidth=1)
ax.text(815, 2.6, 'S2 cleavage site 815-816', ha='center', va='bottom', color='black', fontsize=10)


# Show negative y-axis ticks
ax.set_yticks(range(-5, 6))
#ax.set_yticklabels(range(-5, 6)) # This part is to show the Y axis scale - the 3 lines below remove elemets of the scale one by one
ax.set_yticklabels([])  # Remove y-axis scale (numbers)
ax.spines['left'].set_linewidth(0)  # Remove y-axis line 
ax.tick_params(axis='y', left=False, length=0) # removes the tiny horizontal lines that demacate the scale

# Ensure zero lines of x and y axes intersect
ax.spines['bottom'].set_position('zero')
ax.spines['left'].set_position('zero')

# Move x-axis label to below the end of the maximum negative y-axis tick
ax.xaxis.set_label_coords(0.5, -0.08)
ax.set_xlabel('Amino Acid Position on Spike Gene')

#Remove borders around figure
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)


# Decrease font size of x-axis tick labels
plt.xticks(fontsize=7)

# Extend positive y-axis limit
ax.set_ylim(-8, 8)  # Adjust as needed

# Plot settings
ax.set_xlim(0, 1273)  # Adjust x-axis limit
ax.set_ylabel('')  # Clear y-axis label
# ax.set_title('Amino Acid Mutations on Spike Gene with Regions')
ax.grid(axis='x', linestyle='--', alpha=0)
plt.tight_layout()

# Add legend at the bottom right corner
handles = [Line2D([0], [0], color=color, lw=5) for color in region_colors.values()]
labels = list(region_full_names.values())
ax.legend(handles, labels, loc='lower right', fontsize='medium', frameon=False)

plt.savefig("spike_mut_40%_REV_BA322.svg")
plt.show()