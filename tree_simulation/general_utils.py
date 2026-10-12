import math
from pathlib import Path

def log_and_save_parameters(cfg: dict, output_dir: Path, filename: str = "parameters.txt"):
    """
    Prints configuration parameters to stdout and writes them to a text file.
    
    Args:
        cfg (dict): Configuration dictionary containing parsed settings.
        output_dir (Path): Output directory path where the log file will be saved.
        filename (str): Name of the text file (default: 'parameters.txt').
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    log_file_path = output_dir / filename

    print("\n" + "=" * 50)
    print(f" SIMULATION PARAMETERS ({output_dir.name})")
    print("=" * 50)

    lines = [f"Simulation Parameters - Output Directory: {output_dir}\n", "=" * 50]

    for key, val in cfg.items():
        # Optional: Format special values like math.inf cleanly
        formatted_val = "inf" if val == math.inf else val
        
        # Print to console
        print(f"  {key:<25}: {formatted_val}")
        
        # Append line for file output
        lines.append(f"{key}: {formatted_val}")

    print("=" * 50 + "\n")

    # Write parameters to text file
    with open(log_file_path, "w") as f:
        f.write("\n".join(lines) + "\n")

    print(f"Saved parameter record to: {log_file_path}")