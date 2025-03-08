import os
import sys
import subprocess
import shutil
import json


def main():
    wizard_root = sys.argv[1]
    env_name = sys.argv[2]

    with open(os.path.join(wizard_root, "installation_data.json"), "r") as f:
        data = json.load(f)
        installed_wizard_dir = data["installed_wizard_dir"]

    shutil.copy(
        "installation_data.json",
        os.path.join(installed_wizard_dir, "installation_data.json"),
    )

    installer_path = os.path.join(
        wizard_root,
        "ext",
        "pymol_wizard_installer",
        "pymol_wizard_installer",
        "install_wizard.py",
    )
    paratope_wizard_root = os.path.join(wizard_root, "ext", "paratope_heatmap_wizard")
    try:
        subprocess.run(
            [
                f"conda run --no-capture-output -n {env_name} python3 {installer_path} {paratope_wizard_root}"
            ],
            shell=True,
            check=True,
        )
    except subprocess.CalledProcessError as e:
        print(f"Something went wrong while installing the paratope heatmap wizard: {e}")
        exit(1)


if __name__ == "__main__":
    main()
