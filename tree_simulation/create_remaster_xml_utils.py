'''
Helper functions to construct ReMaster birth-death tree simulation XML files.

Example XML structure (single tree):

<beast ...>
    <run spec="Simulator" nSims="N">
        <simulate id="tree" spec="SimulatedTree">
            <trajectory id="trajectory" spec="StochasticTrajectory", maxTime="T">
                <population       id="X" spec="RealParameter" value="init_pop"/>
                <samplePopulation id="S" spec="RealParameter" value="0"/>
                <samplePopulation id="removed" spec="RealParameter" value="0"/>
                
                <reaction spec="Reaction" rate="lambda"> X -> 2X </reaction>
                <reaction spec="Reaction" rate="mu"> X -> removed </reaction>
                <reaction spec="Reaction" rate="psi"> X -> S + removed </reaction>
            </trajectory>
        </simulate>
        <logger fileName="base.traj"> <log idfer="trajectory"/> </logger>
        <logger fileName="base.trees"> <log idfer="tree"/> </logger>
    </run>         
</beast>
'''

# Indentation helper function
def indent_n_tabs(n=1):
    """
    Build a string of tab characters, used to control indentation depth
    when assembling multi-line XML blocks elsewhere in this module.

    Args:
        n (int, optional): Number of tab characters to produce. Defaults to 1.

    Returns:
        str: A string consisting of `n` tab ("\\t") characters.
    """
    return "\t" * n

# Population initialization helper functions
def initialize_population(name: str, value: int=1) -> str:
    """
    Build a single <population> XML element representing a population
    compartment (e.g. the living/infected population "X") with a given
    starting size.

    Args:
        name (str): The `id` to assign to the population (e.g. "X").
        value (int, optional): Initial population size. Defaults to 1.

    Returns:
        str: A single-line `<population .../>` XML tag.
    """
    return f'<population id="{name}" spec="RealParameter" value="{value}"/>'

def initialize_sample_population(name: str, value: int=0) -> str:
    """
    Build a single <samplePopulation> XML element. Sample populations track
    individuals removed from the simulation via a named "sink" reaction,
    e.g. sampled individuals ("S") or unsampled removals ("removed").

    Args:
        name (str): The `id` to assign to the sample population (e.g. "S").
        value (int, optional): Initial size, almost always 0. Defaults to 0.

    Returns:
        str: A single-line `<samplePopulation .../>` XML tag.
    """
    return f'<samplePopulation id="{name}" spec="RealParameter" value="{value}"/>'


# Reaction helper function
def create_reaction(rate: list[float], left_side: str="X", right_side: str="2X", change_times: list[float] | None=None) -> str:
    """
    Build a single <reaction> XML element describing a stochastic event
    rule (e.g. birth, death, or sampling) of the form
    `left_side -> right_side`, occurring at the given rate. Supports
    piecewise-constant (time-varying) rates via `change_times`.

    Args:
        rate (list[float]): One rate value for a constant rate, or multiple
            values paired with `change_times` for a piecewise rate
            schedule (should have length `len(change_times) + 1`).
        left_side (str, optional): The reactant side of the reaction.
            Defaults to "X".
        right_side (str, optional): The product side of the reaction.
            Defaults to "2X".
        change_times (list[float] | None, optional): Time points at which
            the rate changes. If provided, a `changeTimes` attribute is
            added to the generated tag. Defaults to None.

    Returns:
        str: A single-line `<reaction .../>` XML tag.
    """
    rate_str = " ".join([str(r) for r in rate])
    if (change_times is not None) and (change_times != []):
        change_times_str = " ".join([str(t) for t in change_times])
        return f'<reaction spec="Reaction" rate="{rate_str}" changeTimes="{change_times_str}"> {left_side} -> {right_side} </reaction>'
    else:
        return f'<reaction spec="Reaction" rate="{rate_str}"> {left_side} -> {right_side} </reaction>'

def create_punctual_reaction(probabilities: list[float], times: list[float], left_side: str="X", right_side: str="removed") -> str:
    probabilities_str = " ".join([str(p) for p in probabilities])
    times_str = " ".join([str(t) for t in times])
    
    return f'<reaction spec="PunctualReaction" p="{probabilities_str}" times="{times_str}"> {left_side} -> {right_side} </reaction>'

# Trajectory helper function
def create_trajectory(populations: list, 
                      reactions: list,
                      punctual_reactions: list | None=None,
                      id: str="trajectory", 
                      spec: str="StochasticTrajectory", 
                      max_time: float=10.0,
                      indent: int=0,
                      ends_when: str | None=None,
                      must_have: str | None=None
                      ) -> str:
    """
    Assemble a full <trajectory> block: a list of population / sample
    population elements followed by a list of reaction elements, wrapped
    in a <trajectory> tag with optional stopping conditions.

    Args:
        populations (list): Pre-built `<population>` / `<samplePopulation>`
            XML strings (e.g. produced by `initialize_population` and
            `initialize_sample_population`).
        reactions (list): Pre-built `<reaction>` XML strings (e.g. produced
            by `create_reaction`).
        id (str, optional): The `id` attribute of the `<trajectory>`
            element. Defaults to "trajectory".
        spec (str, optional): The BEAST2/ReMaster `spec` class for the
            trajectory. Defaults to "StochasticTrajectory".
        max_time (float, optional): Maximum simulated time (`maxTime`
            attribute). Defaults to 10.0.
        indent (int, optional): Number of tabs to indent the outer
            `<trajectory>` tag; inner elements are indented one level
            deeper. Defaults to 0.
        ends_when (str | None, optional): Optional expression for the
            `endsWhen` attribute (a simulation stopping condition).
            Defaults to None.
        must_have (str | None, optional): Optional expression for the
            `mustHave` attribute (a condition the simulation must satisfy,
            e.g. "S>=1"). Defaults to None.

    Returns:
        str: A multi-line, indented `<trajectory>...</trajectory>` block.
    """
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
    
    if punctual_reactions is not None:
        for punctual_react in punctual_reactions:
            trajectory += indent_n_tabs(indent+1) + punctual_react + "\n"
            
    trajectory += indent_n_tabs(indent) + '</trajectory>'
    
    
    
    return trajectory

# Simulated tree helper function
def create_tree_simulation(trajectory: str, tree_id: str="tree", indent: int=0) -> str:
    """
    Wrap a pre-built <trajectory> block inside a
    <simulate spec="SimulatedTree"> element.

    Args:
        trajectory (str): The XML string produced by `create_trajectory`.
        tree_id (str, optional): The `id` to assign to the simulated tree.
            Defaults to "tree".
        indent (int, optional): Indentation level for the `<simulate>` tag.
            Defaults to 0.

    Returns:
        str: A multi-line `<simulate>...</simulate>` block containing the
            trajectory.
    """
    simulation = f'{indent_n_tabs(indent)}<simulate id="{tree_id}" spec="SimulatedTree">\n'
    simulation += trajectory + "\n"
    simulation += indent_n_tabs(indent) + '</simulate>\n'

    return simulation

# Logger elements helper function
def create_loggers(basefile_name: str,
                   traj_id: str,
                   tree_id: str,
                   samp_name: str="S",
                   removed_name: str="removed",
                   indent: int=0) -> str:
    """
    Build the three standard <logger> elements used to write simulation
    output to disk: a trajectory logger, a "full tree" logger pruned to
    the removed (unsampled) population, and a "sampled tree" logger
    pruned to the sampled population.

    Args:
        basefile_name (str): Base file name/path; extensions
            (.traj, .full.trees, .sampled.trees) are appended
            automatically.
        traj_id (str): `id` of the trajectory element to log (referenced
            via `idref`).
        tree_id (str): `id` of the simulated tree element to log
            (referenced via `simulatedTree`).
        samp_name (str, optional): Name of the sampled population used for
            the sampled-tree logger. Defaults to "S".
        removed_name (str, optional): Name of the removed population used
            for the full-tree logger. Defaults to "removed".
        indent (int, optional): Indentation level for the logger elements.
            Defaults to 0.

    Returns:
        str: Three newline-separated `<logger>...</logger>` blocks for
            "{basefile_name}.traj", "{basefile_name}.full.trees", and
            "{basefile_name}.sampled.trees".
    """
    logger_trajectory = f'{indent_n_tabs(indent)}<logger fileName="{basefile_name}.traj"> <log idref="{traj_id}"/> </logger>'
    
    logger_full_tree = f'{indent_n_tabs(indent)}<logger fileName="{basefile_name}.full.trees" mode="tree">\n'
    logger_full_tree += f'{indent_n_tabs(indent+1)}<log spec="TypedTreeLogger"> <tree spec="PrunedTree" samplePops="{removed_name}" simulatedTree="@{tree_id}"/> </log>\n'
    logger_full_tree += f'{indent_n_tabs(indent)}</logger>'
    
    logger_sampled_tree = f'{indent_n_tabs(indent)}<logger fileName="{basefile_name}.sampled.trees" mode="tree">\n'
    logger_sampled_tree += f'{indent_n_tabs(indent+1)}<log spec="TypedTreeLogger"> <tree spec="PrunedTree" samplePops="{samp_name}" simulatedTree="@{tree_id}"/> </log>\n'
    logger_sampled_tree += f'{indent_n_tabs(indent)}</logger>'
    
    return logger_trajectory + "\n" + logger_full_tree + "\n" + logger_sampled_tree

# Write full XML file
def write_single_tree_xml(file_path: str, 
                          n_sims: int, 
                          birth_rate: list[float], 
                          death_rate: list[float],
                          change_times_birth: list[float] | None=None,
                          change_times_death: list[float] | None=None,
                          change_times_sampling: list[float] | None=None,
                          max_time: float=10.,
                          sampling_rate: list[float]=[1.], 
                          loggers_base="simulation",
                          traj_type: str="StochasticTrajectory",
                          init_pop: int=1,
                          ends_when: str | None=None,
                          must_have: str | None="S>=1"
                          ) -> None:
    """
    Build a complete, ready-to-run ReMaster BEAST2 XML file for simulating
    a single birth-death tree, and write it to disk.

    Internally this initializes three populations (X: living, S: sampled,
    removed: unsampled removals), builds birth/death/sampling reactions,
    wraps them in a trajectory and tree simulation, attaches the standard
    trajectory/full-tree/sampled-tree loggers, and writes the resulting
    <beast>...</beast> document to `file_path`.

    Args:
        file_path (str): Output path for the generated .xml file. Any
            existing file at this path is overwritten.
        n_sims (int): Number of simulation replicates (`nSims` attribute).
        birth_rate (list[float]): Birth rate(s); a single-element list for
            a constant rate, or multiple values for a piecewise-constant
            rate schedule paired with `change_times_birth`.
        death_rate (list[float]): Death rate(s); same convention as
            `birth_rate`, paired with `change_times_death`.
        change_times_birth (list[float] | None, optional): Time points at
            which the birth rate changes. Defaults to None.
        change_times_death (list[float] | None, optional): Time points at
            which the death rate changes. Defaults to None.
        change_times_sampling (list[float] | None, optional): Time points
            at which the sampling rate changes. Defaults to None.
        max_time (float, optional): Maximum simulated time. Defaults to 10.
        sampling_rate (list[float], optional): Sampling rate(s); same
            convention as `birth_rate`, paired with `change_times_sampling`.
            Defaults to [1.].
        loggers_base (str, optional): Base name used for the three output
            log files. Defaults to "simulation".
        traj_type (str, optional): ReMaster trajectory `spec` class to use.
            Defaults to "StochasticTrajectory".
        init_pop (int, optional): Initial size of population X. Defaults to 1.
        ends_when (str | None, optional): Optional simulation-ending
            condition expression. Defaults to None.
        must_have (str | None, optional): Condition the simulation must
            satisfy; defaults to "S>=1" (at least one sample).

    Returns:
        None: Writes the XML file to `file_path` as a side effect.
    """
    populations = [initialize_population("X", init_pop), 
                   initialize_sample_population("S", 0),
                   initialize_sample_population("removed", 0)]
    reactions = [create_reaction(birth_rate, "X", "2X", change_times_birth), 
                 create_reaction(death_rate, "X", "removed", change_times_death), 
                 create_reaction(sampling_rate, "X", "S + removed", change_times_sampling)]
    punctual_reactions = [create_punctual_reaction([1.0], max_time, "X", "removed")]
    
    with open(file_path, "w") as f:
        f.write('<beast version="2.0" namespace="beast.base.inference:beast.base.inference.parameter:remaster">\n')
        f.write(f'\t<run spec="Simulator" nSims="{n_sims}">\n')
        f.write(create_tree_simulation(create_trajectory(populations, 
                                                         reactions,
                                                         punctual_reactions, 
                                                         max_time=max_time, 
                                                         spec=traj_type, 
                                                         ends_when=ends_when, 
                                                         must_have=must_have, 
                                                         indent=3), indent=2) + "\n")
        f.write(create_loggers(loggers_base, "trajectory", "tree", "S", "removed", indent=2) + "\n")
        f.write('\t</run>\n')
        f.write('</beast>\n')
        
# Write full XML file
def write_single_tree_for_wastewater_xml(file_path: str, 
                                         n_sims: int, 
                                         birth_rate: list[float], 
                                         death_rate: list[float],
                                         change_times_birth: list[float] | None=None,
                                         change_times_death: list[float] | None=None,
                                         change_times_sampling: list[float] | None=None,
                                         max_time: float=10.,
                                         sequencing_delay: float=0.,
                                         sampling_rate: list[float]=[1.], 
                                         loggers_base="simulation",
                                         traj_type: str="StochasticTrajectory",
                                         init_pop: int=1,
                                         ends_when: str | None=None,
                                         must_have: str | None="S>=1"
                                         ) -> None:
    """
    Build a complete, ready-to-run ReMaster BEAST2 XML file for simulating
    a single birth-death tree, and write it to disk.

    Internally this initializes three populations (X: living, S: sampled,
    removed: unsampled removals), builds birth/death/sampling reactions,
    wraps them in a trajectory and tree simulation, attaches the standard
    trajectory/full-tree/sampled-tree loggers, and writes the resulting
    <beast>...</beast> document to `file_path`.

    Args:
        file_path (str): Output path for the generated .xml file. Any
            existing file at this path is overwritten.
        n_sims (int): Number of simulation replicates (`nSims` attribute).
        birth_rate (list[float]): Birth rate(s); a single-element list for
            a constant rate, or multiple values for a piecewise-constant
            rate schedule paired with `change_times_birth`.
        death_rate (list[float]): Death rate(s); same convention as
            `birth_rate`, paired with `change_times_death`.
        change_times_birth (list[float] | None, optional): Time points at
            which the birth rate changes. Defaults to None.
        change_times_death (list[float] | None, optional): Time points at
            which the death rate changes. Defaults to None.
        change_times_sampling (list[float] | None, optional): Time points
            at which the sampling rate changes. Defaults to None.
        max_time (float, optional): Maximum simulated time. Defaults to 10.
        sampling_rate (list[float], optional): Sampling rate(s); same
            convention as `birth_rate`, paired with `change_times_sampling`.
            Defaults to [1.].
        loggers_base (str, optional): Base name used for the three output
            log files. Defaults to "simulation".
        traj_type (str, optional): ReMaster trajectory `spec` class to use.
            Defaults to "StochasticTrajectory".
        init_pop (int, optional): Initial size of population X. Defaults to 1.
        ends_when (str | None, optional): Optional simulation-ending
            condition expression. Defaults to None.
        must_have (str | None, optional): Condition the simulation must
            satisfy; defaults to "S>=1" (at least one sample).

    Returns:
        None: Writes the XML file to `file_path` as a side effect.
    """
    sequencing_delay_in_years = sequencing_delay / 365.25
    end_of_simulation = max_time + sequencing_delay_in_years
    sampling_rate.append(0)
    if change_times_sampling is not None:
        change_times_sampling.append(max_time)
    else:
        change_times_sampling = [max_time]
    
    populations = [initialize_population("X", init_pop), 
                   initialize_sample_population("S", 0),
                   initialize_sample_population("removed", 0)]
    reactions = [create_reaction(birth_rate, "X", "2X", change_times_birth), 
                 create_reaction(death_rate, "X", "removed", change_times_death), 
                 create_reaction(sampling_rate, "X", "S + removed", change_times_sampling)]
    punctual_reactions = [create_punctual_reaction([1.0], [end_of_simulation], "X", "removed")]
    
    with open(file_path, "w") as f:
        f.write('<beast version="2.0" namespace="beast.base.inference:beast.base.inference.parameter:remaster">\n')
        f.write(f'\t<run spec="Simulator" nSims="{n_sims}">\n')
        f.write(create_tree_simulation(create_trajectory(populations, 
                                                         reactions,
                                                         punctual_reactions, 
                                                         max_time=end_of_simulation, 
                                                         spec=traj_type, 
                                                         ends_when=ends_when, 
                                                         must_have=must_have, 
                                                         indent=3), indent=2) + "\n")
        f.write(create_loggers(loggers_base, "trajectory", "tree", "S", "removed", indent=2) + "\n")
        f.write('\t</run>\n')
        f.write('</beast>\n')



# Only for testing purposes
# def main():
#     file_path = "./test.xml"
#     n_sims = 1
#     birth_rate = [1.4, 1.2, 1.0]
#     death_rate = [0.5, 0.4, 0.3]
#     sampling_rate = [0.5, 0.4, 0.3]   
#     change_times = [2.0, 5.7]
#     write_single_tree_xml(file_path, n_sims, birth_rate, death_rate, change_times=change_times, max_time=10.0, sampling_rate=sampling_rate, loggers_base="simulation", traj_type="StochasticTrajectory", init_pop=1)
    
#     # populations = [initialize_population("X"), initialize_sample_population("S")]
#     # reactions = [create_reaction(birth_rate, "X", "2X"), create_reaction(death_rate, "X", "0"), create_reaction(sampling_rate, "X", "S")]
#     # trajectory = create_trajectory(populations, reactions, max_time=10, indent=1)
#     # simulated_tree = create_tree_simulation(trajectory)
#     # loggers = create_loggers("simulation", "trajectory", "tree", "S", indent=0)
    
#     # with open("test.xml", "w") as f:
#     #     f.write(simulated_tree)
#     #     f.write("\n")
#     #     f.write(loggers)

# if __name__ == "__main__":
#     main()
