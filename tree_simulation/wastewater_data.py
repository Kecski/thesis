import numpy as np
import pandas as pd

traj_file = "./remaster_sims/logs/sim.traj"

traj = pd.read_csv(traj_file, sep='\t')

print(traj)