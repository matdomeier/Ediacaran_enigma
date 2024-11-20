## Ediacaran_enigma
This repository includes a set of codes that can be used to generate the simulations shown and described in:
Domeier, M., Robert, B., Meert, J.G., Kulakov, E.V., McCausland, P.J.A., Trindade, R.I.F. and Torsvik, T.H. (2023) The enduring Ediacaran paleomagnetic enigma, Earth-Science Reviews, 242, 104444, DOI:10.1016/j.earscirev.2023.104444.

#### Software environment
The codes were written in Python using Jupyter notebooks. The versions of the packages used are provided in the `environment.yml` file, which can also be used to recreate the original conda environment used.

#### Description of contents

`make_toy_models.ipynb` is the main notebook, which draws information from the other files provided. With this notebook, a base kinematic simulation can be executed, as well as 4 derivative simulations that model: 1) ultra-fast plate drift, 2) data-corruption, 3) true polar wander, and 4) an anomalous magnetic field.

`auxiliary_functions.py` provides a set of functions needed to execute various mathematical operations necessary in `make_toy_models.ipynb`.

`arbitrary_polygons.gpml` provides the set of 3 arbitrary polygons used in the original publication.

Within the rot_models directory are a set of rotation files that provide the basis kinematics for the 4 alternative derivative models. These .rot files have the following column structure: moving_plate_ID, time, euler_pole_lat, euler_pole_lon, euler_pole_angle, reference_plate_id.