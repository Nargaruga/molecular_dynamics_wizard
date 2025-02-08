import setuptools

setuptools.setup(
  name="molecular_dynamics_wizard", 
  version="0.1.0",          
  packages=["molecular_dynamics"],  
  entry_points = {
      'console_scripts': ['mdsim=molecular_dynamics.simulator:main'],
  }
)
