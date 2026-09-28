import numpy as np
import pandas as pd
from datetime import datetime, timedelta
from pathlib import Path

length_in_years = 0.5
length_in_days = length_in_years * 365.25
traj_file = "./remaster_sims/logs/sim.traj"
start_date = datetime(2002, 1, 1)

def get_pop_trajectory(traj_file: str, pop :str="X") -> pd.DataFrame:
    traj = pd.read_csv(traj_file, sep='\t')
    traj_pop = traj.loc[traj.loc[:, "population"] == pop, ["t", "value"]].reset_index(drop=True) # Trajectory of the given population
    
    return traj_pop

def add_noise(measures: np.array, std: float=1, seed: int=42) -> np.array:
    rng = np.random.default_rng(seed=seed)
    noise = rng.normal(loc=0, scale=std, size=measures.shape)
    
    return measures + noise


def create_wastewater_concentrations(traj: pd.DataFrame,
                                  output_dir: str,
                                  sim_length: float=1.0,
                                  start_date: datetime=(2002, 1, 1),
                                  delay: float=0.0,
                                  noise_std: float=1.0
                                  ) -> None:
    '''
    sim_length (years)
    delay (days)
    '''
    sim_length_in_days = sim_length * 365.25
    delay_in_years = delay / 365.25
    
    measure_times = np.arange(0, sim_length, sim_length/sim_length_in_days)
    measure_datetimes = [start_date + timedelta(days=years * 365.25) for years in measure_times] 
    measure_dates = [d.date() for d in measure_datetimes]
    
    measures = []
    for t in measure_times:
        idx = traj["t"].searchsorted(t-delay_in_years, side="left") - 1
        if idx >= 0:
            measures.append(traj.loc[idx, "value"])
        else:
            measures.append(0)
            
    noisy_measures = add_noise(np.array(measures), std=noise_std)
    
    measures_df = pd.DataFrame({"date": measure_dates,
                                    "concentration": noisy_measures})
    flows_df = pd.DataFrame({"date": measure_dates,
                                "flow": 1})
    
    output_dir_path = Path(output_dir)
    output_dir_path.mkdir(parents=True, exist_ok=True)
    
    measures_df.to_csv(f"{output_dir}/concentrations.csv", index=False)
    flows_df.to_csv(f"{output_dir}/flows.csv", index=False)
    

def main():
    length_in_years = 0.5
    # length_in_days = length_in_years * 365.25
    traj_file = "./remaster_sims/logs/sim.traj"
    output_dir = "./episewer_learn"
    start_date = datetime(2002, 1, 1)
    delay_in_days = 7.0
    
    traj_X = get_pop_trajectory(traj_file=traj_file, pop="X")
    create_wastewater_concentrations(traj_X, output_dir, length_in_years, start_date, delay=delay_in_days)
    
    
    # measure_times = np.arange(0, length_in_years, length_in_years/length_in_days) # So we have wastewater measurements every day

    # # traj = pd.read_csv(traj_file, sep='\t')
    # # traj_X = traj.loc[traj.loc[:, "population"] == "X", ["t", "value"]].reset_index(drop=True) # Trajectory of the whole population

    # # Creating wastewater measurements with added random noise
    # measures = []
    # for t in measure_times:
    #     idx = traj_X["t"].searchsorted(t, side="left") - 1
    #     if idx >= 0:
    #         measures.append(traj_X.loc[idx, "value"])
    #     else:
    #         measures.append(0)
        
    # measures = np.array(measures)

    # # Adding noise
    # rng = np.random.default_rng(seed=42)
    # noise = rng.normal(loc=0, scale=0.5, size=measures.shape)
    # noisy_measures = measures + noise

    # measure_datetimes = [start_date + timedelta(days=years * 365.25) for years in measure_times] 
    # measure_dates = [d.date() for d in measure_datetimes]

    # print(f"measures: {measures}")
    # print(f"noisy measures: {noisy_measures}")

    # measures_df = pd.DataFrame({"date": measure_dates,
    #                             "concentration": noisy_measures})
    # flows_df = pd.DataFrame({"date": measure_dates,
    #                         "flow": 1})

    # measures_df.to_csv("./episewer_learn/concentrations.csv", index=False)
    # flows_df.to_csv("./episewer_learn/flows.csv", index=False)
    

if __name__ == "__main__":
    main()