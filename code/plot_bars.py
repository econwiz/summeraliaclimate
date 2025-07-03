from pathlib import Path
import matplotlib.pyplot as plt

gwl = [1.0, 1.5, 2.0, 2.5, 3.0]

f_internal_mort = [0.575697, 0.440232, 0.208426, 0.161265, 0.107334]
f_model_mort    = [0.424303, 0.559768, 0.791574, 0.838735, 0.892666]

f_internal_gdp  = [0.502662, 0.535638, 0.408623, 0.316119, 0.218661]
f_model_gdp     = [0.497338, 0.464362, 0.591377, 0.683881, 0.781339]

project_root = Path(__file__).resolve().parent.parent
figs_dir     = project_root / "figs"
figs_dir.mkdir(parents=True, exist_ok=True)

fig1, ax1 = plt.subplots()
ax1.bar(gwl, f_model_mort,    width=0.4, label='model')
ax1.bar(gwl, f_internal_mort, width=0.4, bottom=f_model_mort, label='internal')
ax1.set(xlabel='Global-warming level (°C)', ylabel='Fraction of variance',
        title='Mortality uncertainty fractions by GWL')
ax1.set_xticks(gwl)
ax1.legend()
plt.tight_layout()
fig1.savefig(figs_dir / 'mortality_bar_by_gwl.png', dpi=300, bbox_inches='tight')
plt.close(fig1)

# 5) Plot & save GDP
fig2, ax2 = plt.subplots()
ax2.bar(gwl, f_model_gdp,    width=0.4, label='model')
ax2.bar(gwl, f_internal_gdp, width=0.4, bottom=f_model_gdp, label='internal')
ax2.set(xlabel='Global-warming level (°C)', ylabel='Fraction of variance',
        title='Per-capita GDP uncertainty fractions by GWL')
ax2.set_xticks(gwl)
ax2.legend()
plt.tight_layout()
fig2.savefig(figs_dir / 'gdp_bar_by_gwl.png', dpi=300, bbox_inches='tight')
plt.close(fig2)