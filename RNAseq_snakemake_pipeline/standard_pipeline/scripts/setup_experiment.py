import yaml
import subprocess

config_file = "RNAseq.standard.Snakemake.config.yaml"

# Load existing config
with open(config_file, "r") as f:
    config = yaml.safe_load(f)

# Ask user if they want to run DESeq2
run_deseq2 = input("Do you want to run DESeq2 analysis? (yes/no): ").strip().lower()

if run_deseq2 == "yes":
    control_samples = input("Enter control sample names (comma-separated): ").strip().split(",")
    treatment_samples = input("Enter treatment sample names (comma-separated): ").strip().split(",")

    config["deseq2"] = {
        "enabled": True,
        "control": control_samples,
        "treatment": treatment_samples
    }
else:
    config["deseq2"] = {"enabled": False}

# Save updated config
with open(config_file, "w") as f:
    yaml.dump(config, f)

print("Config updated! Running Snakemake...")

# Run Snakemake
try:
    subprocess.run(["snakemake", "--cores", "4"], check=True)
except subprocess.CalledProcessError as e:
    print(f"Error running Snakemake: {e}")