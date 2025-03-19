import os
import json
from pathlib import Path

from pymol.Qt import QtWidgets

from .settings_gui import Ui_Form
from molecular_dynamics.simulation_params import SimulationParameters

plugin_root = Path(__file__).parent

# global reference to avoid garbage collection of our dialog
dialog = None


def run_plugin_gui():
    global dialog

    if dialog is None:
        dialog = QtWidgets.QDialog()
        form = Ui_Form()
        form.setupUi(dialog)
        load_current_parameters(form)

        form.apply_btn.clicked.connect(lambda: save_parameters(form))

    dialog.show()


def get_configuration_path():
    install_data_path = os.path.join(plugin_root, "installation_data.json")
    with open(install_data_path) as f:
        data = json.load(f)
        try:
            installed_wizard_dir = data["installed_wizard_dir"]
            return os.path.join(
                installed_wizard_dir, "dynamics_extra", "simulation_params.yaml"
            )
        except KeyError:
            print(
                f"WARNING: Failed to read configuration file path from {install_data_path}."
            )

    return None


def load_configuration():
    params = SimulationParameters()

    config_path = get_configuration_path()
    if config_path is None:
        return params
    else:
        params.parse_file(config_path)

    return params


def save_configuration(params):
    config_path = get_configuration_path()
    if config_path is None:
        return
    else:
        params.serialize(config_path)


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
