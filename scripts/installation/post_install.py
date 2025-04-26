import os
import sys
import subprocess


def main():
    wizard_root = sys.argv[1]
    env_name = sys.argv[2]

    paratope_wizard_root = os.path.join(wizard_root, "ext", "paratope_heatmap_wizard")

    if os.name == "nt":
        prefix = ["powershell.exe"]
    else:
        prefix = []

    try:
        subprocess.run(
            prefix
            + [
                "conda",
                "run",
                "--name",
                "base",
                "--no-capture-output",
                "install_wizard",
                paratope_wizard_root,
                env_name,
            ],
            check=True,
        )
    except subprocess.CalledProcessError as e:
        print(f"Something went wrong while installing the paratope heatmap wizard: {e}")
        exit(1)

    if os.name == "nt":
        subprocess.run(
            [
                "powershell.exe",
                "Compress-Archive",
                "-Path",
                "wizard_settings_plugin",
                "-DestinationPath",
                "plugin.zip",
                "-Force",
            ],
            cwd=wizard_root,
            check=True,
        )
    else:
        subprocess.run(
            [
                "zip",
                "-r",
                "plugin.zip",
                "wizard_settings_plugin",
            ],
            cwd=wizard_root,
            check=True,
        )


if __name__ == "__main__":
    main()
