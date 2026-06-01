"""Plot the LDA and Koopmans band structures of silicon."""

import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np

from koopmans.io import read

# Load the workflow
wf = read('si.pkl')

# Extract the band structures from the workflow
# Note: the API is in the process of being refactored, hence the difference
# in how we access the band structures for the LDA and Koopmans calculations
koopmans_valence_bs = wf.processes[-2].outputs.band_structure
koopmans_conduction_bs = wf.processes[-1].outputs.band_structure
lda_bs = wf.calculations[-2].results['band structure']

reference = lda_bs.energies[:, :, :4].max()

# Fetch the LDA bands, and shift them by the same amount
lda_bs_shifted = lda_bs.subtract_reference(reference)
koopmans_valence_bs = koopmans_valence_bs.subtract_reference(reference)
koopmans_conduction_bs = koopmans_conduction_bs.subtract_reference(reference)

def annotate_gap(ax, x_v, y_v, x_c, y_c):
    """Draw a double-arrow between (x_v, y_v) and (x_c, y_c) labelled with the gap in eV."""
    ax.annotate(xy=(x_v, y_v), xycoords='data',
                xytext=(x_c, y_c), textcoords='data', text='',
                arrowprops={'arrowstyle': '<->', 'shrinkA': 0, 'shrinkB': 0})
    ax.text((x_c + x_v) / 2, (y_c + y_v) / 2,
            f'{y_c - y_v:.2f} eV', ha='center', va='center',
            path_effects=[pe.withStroke(linewidth=4, foreground='white')])


# Plot the two band structures
ax = lda_bs_shifted.plot(label='LDA', spin=0, color='tab:blue', ls='--')
ax = koopmans_valence_bs.plot(ax=ax, label='KI@LDA', color='tab:green')
ax = koopmans_conduction_bs.plot(ax=ax, color='tab:green')

x, label_xcoords, labels = koopmans_valence_bs.get_labels()

# Find the Koopmans valence band maximum and conduction band minimum
valence = koopmans_valence_bs.energies
conduction = koopmans_conduction_bs.energies
i_vbm = np.unravel_index(np.nanargmax(valence), valence.shape)
i_cbm = np.unravel_index(np.nanargmin(conduction), conduction.shape)

# Label the fundamental band gap
annotate_gap(ax, x[i_vbm[1]], valence[i_vbm], x[i_cbm[1]], conduction[i_cbm])

# Find and label the Koopmans Gamma -> Gamma direct gap
i_gamma = next(i for i, lbl in enumerate(labels) if lbl in ('G', 'Г', r'$\Gamma$', 'Gamma'))
x_gamma = label_xcoords[i_gamma]
i_kpt_gamma = int(np.argmin(np.abs(x - x_gamma)))
y_gamma_v = np.nanmax(valence[:, i_kpt_gamma, :])
y_gamma_c = np.nanmin(conduction[:, i_kpt_gamma, :])
annotate_gap(ax, x_gamma, y_gamma_v, x_gamma, y_gamma_c)
 
# Tweak the figure aesthetics
ax.legend(loc='lower right', ncol=2, bbox_to_anchor=(1, 1))
ax.set_ylim([-10, 5])

# Display or save the figure (uncomment as desired)
plt.savefig('si_bandstructures.png')
# plt.show()
