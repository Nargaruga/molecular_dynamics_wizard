import os
import shutil
import json


def main():
    with open("installation_data.json", "r") as f:
        data = json.load(f)
        installed_wizard_dir = data["installed_wizard_dir"]

    shutil.copy(
        "installation_data.json",
        os.path.join(installed_wizard_dir, "installation_data.json"),
    )


if __name__ == "__main__":
    main()
