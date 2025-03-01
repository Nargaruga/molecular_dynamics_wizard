import os
from pathlib import Path
import subprocess
import shutil
import re
from util.constants import PYMOL_PYTHON_VERSION, ENTRY_NAME, WIZARD_NAME
import json

wizard_root = Path(__file__).parent.parent


def add_line_after(file, to_insert, pattern_to_insert, target_pattern):
    """Adds a line to the file after the specified point. If the line is already present, the file is unchanged."""

    with open(file, "r") as f:
        contents = f.read()

    if re.search(pattern_to_insert, contents):
        print(f"Entry already exists in {file}, skipping...")
        return

    target = target_pattern.search(contents)
    if target is None:
        print(f"Could not find target in {file}")
        return

    target_end = target.end()

    contents = contents[:target_end] + f"{to_insert}" + contents[target_end:]

    with open(file, "w") as f:
        f.write(contents)


# Entry point
current_env = os.environ.get("CONDA_DEFAULT_ENV")
if current_env is None:
    print("Could not detect conda environment. Is conda installed?")
    exit(1)

print(
    f"You are currently about to install the wizard in the {current_env} environment. Do you wish to continue? (y/N)"
)
try:
    answer = input().strip().lower() or "n"
except KeyboardInterrupt:
    print("Aborted by user.")
    exit(0)

if answer != "y":
    print("Aborted by user.")
    exit(0)

try:
    conda_base_path = str(
        subprocess.check_output("conda info --base", shell=True), "utf-8"
    ).strip()
except subprocess.CalledProcessError:
    print("Failed to retrieve conda base path.")
    exit(1)

prefix = os.path.join(conda_base_path, "envs", current_env)
if prefix is None:
    print("Something went wrong while creating the new environment.")
    exit(1)

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

if not os.path.exists(pymol_dir):
    print(f"PyMOL is not installed in the {current_env} environment. Aborting.")
    exit(1)

print("Copying files...")
installed_wizard_dir = os.path.join(pymol_dir, "wizard")
config_location = os.path.join(prefix, ".pymol_wizards_config", f"{WIZARD_NAME}")

# Record information for the uninstallation process
data = {
    "conda_env": current_env,
    "config_path": config_location,
    "to_remove": [
        os.path.join(installed_wizard_dir, f"{WIZARD_NAME}.py"),
        os.path.join(installed_wizard_dir, f"{WIZARD_NAME}_extra"),
        os.path.join(installed_wizard_dir, "molecular_dynamics"),
    ],
}
with open(os.path.join(wizard_root, "installation_data.json"), "w") as f:
    json.dump(data, f)

try:
    shutil.copy(
        os.path.join(wizard_root, f"{WIZARD_NAME}.py"),
        os.path.join(installed_wizard_dir, f"{WIZARD_NAME}.py"),
    )

    os.makedirs(
        os.path.join(installed_wizard_dir, f"{WIZARD_NAME}_extra"), exist_ok=True
    )

    os.makedirs(os.path.join(config_location), exist_ok=True)

    # shutil.copytree(
    #     os.path.join(wizard_root, "martini_parameters"),

    #     os.path.join(
    #         installed_wizard_dir, f"{WIZARD_NAME}_extra", "martini_parameters"
    #     ),
    #     dirs_exist_ok=True,
    # )

    shutil.copy(
        os.path.join(wizard_root, "simulation_params.yaml"),
        os.path.join(config_location, "simulation_params.yaml"),
    )

    shutil.copy(
        os.path.join(wizard_root, "installation_data.json"),
        os.path.join(
            installed_wizard_dir, f"{WIZARD_NAME}_extra", "installation_data.json"
        ),
    )
    shutil.copy(
        os.path.join(wizard_root, "installation_data.json"),
        os.path.join(wizard_root, "plugin", "installation_data.json"),
    )

    shutil.copytree(
        os.path.join(wizard_root, "molecular_dynamics"),
        os.path.join(installed_wizard_dir, "molecular_dynamics"),
        dirs_exist_ok=True,
    )
except shutil.Error as e:
    print(f"Failed to copy files: {e}")
    exit(1)

settings_plugin_archive = "settings_plugin"
shutil.make_archive(
    settings_plugin_archive, "zip", os.path.join(f"{wizard_root}", "plugin")
)

print("Adding menu entries...")
# Edit the openvr wizard to add a menu item in the internal menu
openvr_wizard = os.path.join(installed_wizard_dir, "openvr.py")
openvr_entry = f"\n[1, '{ENTRY_NAME}', 'wizard {WIZARD_NAME}'],"
openvr_entry_pattern = re.compile(openvr_entry.replace("[", r"\[").replace("]", r"\]"))
openvr_target_pattern = re.compile(r"\[2, 'Wizard Menu', ''\],")
add_line_after(openvr_wizard, openvr_entry, openvr_entry_pattern, openvr_target_pattern)

# Edit the openvr wizard to add a menu item in the internal menu
gui_file = os.path.join(pymol_dir, "_gui.py")
external_entry = f"\n('command', '{ENTRY_NAME}', 'wizard {WIZARD_NAME}'),"
external_entry_pattern = re.compile(
    external_entry.replace("(", r"\(").replace(")", r"\)")
)
external_target_pattern = re.compile(r"\('menu', 'Wizard', \[")
add_line_after(
    gui_file, external_entry, external_entry_pattern, external_target_pattern
)

print("Done!")
print(f"Remember to activate the {current_env} conda environment before running PyMOL.")
print(
    f"The plugin for changing wizard settings is located at {os.path.join(wizard_root, settings_plugin_archive)} and can be installed through PyMOL's plugin manager."
)
