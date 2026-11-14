"""Rerender the three figures from actual compact data; no dynamics."""
from pathlib import Path
import shutil
from coherence_figures import correlation_svg,phase_picture
p=Path(__file__).resolve().parent
cold=['cold_x1_h008_t40', 'cold_x1_h008_L48_t40', 'cold_x1_h008_L100_t40', 'cold_x1_h008_L48_t80']
window=['cold_x1_h008_t40', 'hot_x1_h000_coldclock_t40', 'hot_x1_h008_coldclock_t40', 'hot_x1_h008_raisedclock_t40']
correlation_svg(p/'phase_data/cold',p/'cold_size_time.svg',cold)
correlation_svg(p/'phase_data/windows',p/'temperature_correlations.svg',window,True,True)
phase_picture(p/'phase_data/cold','cold_x1_h008_L100_t40',3.732,True)
shutil.copyfile(p/'phase_data/cold/cold_x1_h008_L100_t40_phase.png',p/'cold_x1_h008_L100_t40_phase.png')
