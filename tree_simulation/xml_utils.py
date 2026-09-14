'''
Helper functions to construct ReMaster birth-death tree simulation XML files.

Example XML structure (single tree):

<beast ...>
    <run spec="Simulator" nSims="N">
        <simulate id="tree" spec="SimulatedTree">
            <trajectory id="trajectory" spec="StochasticTrajectory", maxTime="T">
                <population       id="X" spec="RealParameter" value="init_pop"/>
                <samplePopulation id="X" spec="RealParameter" value="0"/>
                <reaction spec="Reaction" rate="lambda"> X -> 2X </reaction>
                <reaction spec="Reaction" rate="mu"> X -> 0 </reaction>
                <reaction spec="Reaction" rate="psi"> X -> S </reaction>
            </trajectory>
        </simulate>
        <logger fileName="base.traj"> <log idfer="trajectory"/> </logger>
        <logger fileName="base.trees"> <log idfer="tree"/> </logger>
    </run>         
</beast>
'''

# Indentation helper function
def indent_n_tabs(n=1):
    return "\t" * n

# Population initialization helper functions
def initialize_population(name: str, value: int=1) -> str:
    return f'<population id="{name}" spec="RealParameter" value="{value}"/>'

def initialize_sample_population(name: str, value: int=0) -> str:
    return f'<samplePopulation id="{name}" spec="RealParameter" value="{value}"/>'


# Reaction helper function
def create_reaction(rate: float, left_side: str="X", right_side: str="2X") -> str:
    return f'<reaction spec="Reaction" rate="{rate}"> {left_side} -> {right_side} </reaction>'

# Trajectory helper function
def create_trajectory(populations: list, 
                      reactions: list,
                      id: str="trajectory", 
                      spec: str="StochasticTrajectory", 
                      max_time: float=10.0,
                      indent: int=0,
                      ends_when: str | None=None,
                      must_have: str | None=None
                      ) -> str:
    trajectory = f'{indent_n_tabs(indent)}<trajectory id="{id}" spec="{spec}" maxTime="{max_time}"'
    if ends_when is not None:
        trajectory += f' endsWhen="{ends_when}"'
    if must_have is not None:
        trajectory += f' mustHave="{must_have}"'
    trajectory += '>\n'
    
    for pop in populations:
        trajectory += indent_n_tabs(indent+1) + pop + "\n"
    trajectory += "\n"
    
    for react in reactions:
        trajectory += indent_n_tabs(indent+1) + react + "\n"
    trajectory += indent_n_tabs(indent) + '</trajectory>'
    
    return trajectory

# Simulated tree helper function
def create_tree_simulation(trajectory: str, tree_id: str="tree", indent: int=0) -> str:
    simulation = f'{indent_n_tabs(indent)}<simulate id="{tree_id}" spec="SimulatedTree">\n'
    simulation += trajectory + "\n"
    simulation += indent_n_tabs(indent) + '</simulate>\n'

    return simulation

# Logger elements helper function
def create_loggers(basefile_name: str, traj_id: str, tree_id: str, samp_name: str="S", indent: int=0) -> str:
    logger_trajectory = f'{indent_n_tabs(indent)}<logger fileName="{basefile_name}.traj"> <log idref="{traj_id}"/> </logger>'
    logger_full_tree = f'{indent_n_tabs(indent)}<logger fileName="{basefile_name}.full.trees"> <log idref="{tree_id}"/> </logger>'
    logger_sampled_tree = f'{indent_n_tabs(indent)}<logger fileName="{basefile_name}.sampled.trees" mode="tree">\n'
    logger_sampled_tree += f'{indent_n_tabs(indent+1)}<log spec="TypedTreeLogger"> <tree spec="PrunedTree" samplePops="{samp_name}" simulatedTree="@{tree_id}"/> </log>\n'
    logger_sampled_tree += f'{indent_n_tabs(indent)}</logger>'
    
    return logger_trajectory + "\n" + logger_full_tree + "\n" + logger_sampled_tree

# Write full XML file
def write_single_tree_xml(file_path: str, 
                          n_sims: int, 
                          birth_rate: float, 
                          death_rate: float,
                          max_time: float=10.0,
                          sampling_rate: float=1, 
                          loggers_base="simulation",
                          traj_type: str="StochasticTrajectory",
                          init_pop: int=1,
                          ends_when: str | None=None,
                          must_have: str | None=None
                          ) -> None:
    populations = [initialize_population("X", init_pop), initialize_sample_population("S", 0)]
    reactions = [create_reaction(birth_rate, "X", "2X"), create_reaction(death_rate, "X", "0"), create_reaction(sampling_rate, "X", "S")]
    
    with open(file_path, "w") as f:
        f.write('<beast version="2.0" namespace="beast.base.inference:beast.base.inference.parameter:remaster">\n')
        f.write(f'\t<run spec="Simulator" nSims="{n_sims}">\n')
        f.write(create_tree_simulation(create_trajectory(populations, reactions, max_time=max_time, spec=traj_type, ends_when=ends_when, must_have=must_have, indent=3), indent=2) + "\n")
        f.write(create_loggers(loggers_base, "trajectory", "tree", "S", indent=2) + "\n")
        f.write('\t</run>\n')
        f.write('</beast>\n')



# Only for testing purposes
def main():
    file_path = "./test.xml"
    n_sims = 1
    birth_rate = 1.4
    death_rate = 0.5
    sample_rate = 0.5
    write_single_tree_xml(file_path, n_sims, birth_rate, death_rate, max_time=10.0, sampling_rate=sample_rate, loggers_base="simulation", traj_type="StochasticTrajectory", init_pop=1)
    
    # populations = [initialize_population("X"), initialize_sample_population("S")]
    # reactions = [create_reaction(birth_rate, "X", "2X"), create_reaction(death_rate, "X", "0"), create_reaction(sample_rate, "X", "S")]
    # trajectory = create_trajectory(populations, reactions, max_time=10, indent=1)
    # simulated_tree = create_tree_simulation(trajectory)
    # loggers = create_loggers("simulation", "trajectory", "tree", "S", indent=0)
    
    # with open("test.xml", "w") as f:
    #     f.write(simulated_tree)
    #     f.write("\n")
    #     f.write(loggers)

if __name__ == "__main__":
    main()
