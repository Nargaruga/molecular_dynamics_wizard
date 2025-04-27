import os
from pathlib import Path
import pkgutil

from pymol.Qt import QtWidgets

from .settings_gui import Ui_Form
import molecular_dynamics
from molecular_dynamics.simulation_params import SimulationParameters

plugin_root = Path(__file__).parent

# global reference to avoid garbage collection of our dialog
dialog = None

CONFIG_PATH = os.path.join("config", "simulation_params.yaml")


class PluginError(Exception):
    """Custom exception for configuration errors."""

    pass


def run_plugin_gui():
    global dialog

    if dialog is None:
        dialog = QtWidgets.QDialog()
        form = Ui_Form()
        form.setupUi(dialog)
        load_current_parameters(form)

        form.apply_btn.clicked.connect(lambda: save_parameters(form))

    dialog.show()


def load_configuration():
    params = SimulationParameters()

    raw_yaml = pkgutil.get_data("molecular_dynamics", CONFIG_PATH)
    if raw_yaml is None:
        raise PluginError("Could not load simulation parameters.")

    yaml_str = raw_yaml.decode("utf-8")
    params.parse_yaml(yaml_str)

    file_path = os.path.join(Path(molecular_dynamics.__file__).parent, CONFIG_PATH)
    print(f"Loaded configuration from {file_path}")

    return params


def save_configuration(params):
    file_path = os.path.join(Path(molecular_dynamics.__file__).parent, CONFIG_PATH)
    params.serialize(file_path)
    print(f"Saved configuration to {file_path}")


def get_force_field_choices():
    return ["charmm36.xml", "amber14-all.xml"]


def get_water_model_choices():
    return ["charmm36/water.xml", "amber14/tip3p.xml"]


def load_current_parameters(form: Ui_Form):
    params = load_configuration()

    form.add_solvent_checkbox.setChecked(params.add_solvent)
    form.remove_non_sim_checkbox.setChecked(params.remove_non_simulated)

    form.minimization_steps_spin.setValue(params.minimization_steps)

    form.nvt_steps_spin.setValue(params.nvt_steps)
    form.npt_steps_spin.setValue(params.npt_steps)
    form.default_pressure_spin.setValue(params.eq_pressure)
    form.default_temp_spin.setValue(params.eq_temperature)

    form.force_field_combo.addItems(get_force_field_choices())
    form.force_field_combo.setCurrentText(params.force_field)
    form.water_model_combo.addItems(get_water_model_choices())
    form.water_model_combo.setCurrentText(params.water_model)
    form.simulation_steps_spin.setValue(params.sim_steps)

    form.temp_spin.setValue(params.temperature)
    form.friction_coeff_spin.setValue(params.friction_coeff)
    form.time_step_spin.setValue(params.timestep)

    form.report_interval_spin.setValue(params.report_interval)


def save_parameters(form: Ui_Form):
    params = SimulationParameters()

    params.add_solvent = form.add_solvent_checkbox.isChecked()
    params.remove_non_simulated = form.remove_non_sim_checkbox.isChecked()

    params.minimization_steps = form.minimization_steps_spin.value()

    params.nvt_steps = form.nvt_steps_spin.value()
    params.npt_steps = form.npt_steps_spin.value()
    params.eq_pressure = form.default_pressure_spin.value()
    params.eq_temperature = form.default_temp_spin.value()

    params.force_field = form.force_field_combo.currentText()
    params.water_model = form.water_model_combo.currentText()
    params.sim_steps = form.simulation_steps_spin.value()

    params.temperature = form.temp_spin.value()
    params.friction_coeff = form.friction_coeff_spin.value()
    params.timestep = form.time_step_spin.value()

    params.report_interval = form.report_interval_spin.value()

    save_configuration(params)
