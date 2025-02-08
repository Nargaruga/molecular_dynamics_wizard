import os
from pathlib import Path
import subprocess
import shutil
import re
import fileinput
import json
from util.constants import (
    PYMOL_PYTHON_VERSION,
    ENTRY_NAME,
    WIZARD_NAME,
    DEFAULT_ENV_NAME,
)


wizard_root = Path(__file__).parent.parent


def remove_line(file, pattern_to_remove):
    with fileinput.FileInput(file, inplace=True, backup=".bak") as file:
        for line in file:
            if not pattern_to_remove.search(line):
                print(line, end="")


with open(os.path.join(wizard_root, "installation_data.json"), "r") as f:
    data = json.load(f)
    # Retrieve the environment used in the installation
    try:
        env_name = data["conda_env"]
    except KeyError:
        print(
            f'The conda environment used in the installation was not recorded. Please enter the name of the environment, or leave empty for default ("{DEFAULT_ENV_NAME}"):'
        )
        try:
            env_name = input().strip()
            if not env_name:
                env_name = DEFAULT_ENV_NAME
        except KeyboardInterrupt:
            print("Aborted by user.")
            exit(0)

    # Get the files to remove
    try:
        to_remove = data["to_remove"]
    except KeyError:
        print(
            "WARNING: no files recorded for removal, you may have to remove them manually."
        )
        to_remove = []


try:
    conda_base_path = str(
        subprocess.check_output("conda info --base", shell=True), "utf-8"
    ).strip()
except subprocess.CalledProcessError:
    print("Failed to retrieve conda base path.")
    exit(1)

prefix = os.path.join(conda_base_path, "envs", env_name)
if prefix is None:
    print("Something went wrong. Please check the conda environment name.")
    exit(1)

print("Removing files...")

for file in to_remove:
    if os.path.exists(file):
        print(f"Removing {file}...")
        if os.path.isdir(file):
            shutil.rmtree(file)
        else:
            os.remove(file)
    else:
        print(f"File {file} does not exist, skipping...")

print("Removing menu entries...")
if os.name == "nt":
    pymol_dir = os.path.join(
        prefix,
        "Lib",
        "site-packages",
        "pymol",
    )
else:
    pymol_dir = os.path.join(
        prefix,
        "lib",
        f"python{PYMOL_PYTHON_VERSION}",
        "site-packages",
        "pymol",
    )
openvr_wizard_file = os.path.join(pymol_dir, "wizard", "openvr.py")
openvr_entry = f"[1, '{ENTRY_NAME}', 'wizard {WIZARD_NAME}'],\n"
openvr_entry_pattern = re.compile(openvr_entry.replace("[", r"\[").replace("]", r"\]"))
remove_line(openvr_wizard_file, openvr_entry_pattern)

gui_file = os.path.join(pymol_dir, "_gui.py")
external_entry = f"('command', '{ENTRY_NAME}', 'wizard {WIZARD_NAME}'),\n"
external_entry_pattern = re.compile(
    external_entry.replace("(", r"\(").replace(")", r"\)")
)
remove_line(gui_file, external_entry_pattern)

print(
    f"Done! Note that the conda environment will have to be deleted manually. You can do so with the command 'conda env remove -n {env_name}'"
)
