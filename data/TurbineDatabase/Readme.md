# European offshore wind turbine database

A European offshore wind turbine database for mesoscale modelling of European offshore wind farms is presented, containing all operating wind farms as of commissioning date before August 2025. The approach integrates turbine locations from the OpenStreetMap [1] and EMODnet [2] with metadata on turbine and wind farm properties from additional public sources.

The file is called eww_opendatabase.csv. It contains the wind farm name, the Original Equipment Manufacturer, the position (latitude, longitude, country), rated power, rotor diameter, hub height, turbine type and commissioning date of each turbine within the domain.

The turbine information is augmented with generic thrust and power curves calculated via the pyWake [3] turbine generator [4],[5], which are located within the zip-file power_curves.zip.


## Power curves matching the turbine_types in the database csv file

### Infos
To automatically match the turbine_type in the database and the power curve file, ' ' and '/' in the turbine_types column have to be replaced by '_'

Files were created using [4] with the parameters (and default values for all other parameters):

* ws_cutin  = 4.0
* ws_cutout = 25.0
* turbulence_intensity = 0.05

## ​References​
[1] OpenStreetMap contributors (2025). Distributed under the Open Database License (ODbL)

[2] https://emodnet.ec.europa.eu/geoviewer/ ​

[3] Mads M. Pedersen, Alexander Meyer Forsting, Paul van der Laan, Riccardo Riva, Leonardo A. Alcayaga Romàn, Javier Criado Risco, Mikkel Friis-Møller, Julian Quick, Jens Peter Schøler Christiansen, Rafael Valotta Rodrigues, Bjarke Tobias Olsen and Pierre-Elouan Réthoré. (2023, February). PyWake 2.5.0: An open-source wind farm simulation tool. https://gitlab.windenergy.dtu.dk/TOPFARM/PyWake, DTU Wind, Technical University of Denmark.

[4] https://topfarm.pages.windenergy.dtu.dk/PyWake/notebooks/WindTurbines.html ​

[5] https://gitlab.windenergy.dtu.dk/TOPFARM/PyWake/-/blob/fd7a9eb1f6b6d23d07a98717c7f6c106e47fe6f2/py_wake/utils/generic_power_ct_curves.py

See also database_ATTRIBUTION.md

## License
See database_LICENSE_ODbL.txt

## Acknowledgements
The EuroWindWakes European offshore Dataset is developed and released within the EuroWindWakes project using funding proved by CETPartnership, the Clean Energy Transition Partnership under the 2023 joint call for research proposals, co funded by the European Commission (GA 101 069750) and with the funding organizations detailed on https://cetpartnership.eu/funding-agencies-and-call-modules, and co funded by EUDP – The Energy Technology Development and Demonstration Programme, Project nr.: 640245-522075.
