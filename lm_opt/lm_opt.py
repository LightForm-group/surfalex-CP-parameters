'Generate simulations and add their results to an LMFitter object.'

import sys
from pathlib import Path

from ruamel.yaml import YAML

from damask_parse import read_table
from tensile_test import TensileTest, LMFitter, FittingParameter
from tensile_test.utils import read_non_uniform_csv

LM_FIT_JSON_NAME = 'lm_fitter.json'


def read_tensile_test_data(path):
    'Read tensile test data CSV into a list of dict.'

    path = Path(path)
    delimiter = ','

    # Read data:
    headers, arr = read_non_uniform_csv(
        path, delimiter=delimiter, skip_rows=2, header_row=1)

    # Get orientations for each stress/strain column pair:
    with path.open() as handle:
        ln = handle.readline()
        oris = [int(i.split()[0]) for i in ln.split(delimiter)[::2]]

    map_keys = {
        'eng. strain': 'eng_strain',
        'eng. stress /Mpa': 'eng_stress',
    }
    out = [
        {
            'orientation': ori,
            map_keys[headers[ori_idx]]: arr[:, ori_idx],
            map_keys[headers[ori_idx + 1]]: arr[:, ori_idx + 1],
        }
        for ori_idx, ori in enumerate(oris)
    ]

    return out


def set_up():
    'Set up the optimisation process'

    print('{}.set_up'.format(__name__))

    # Load the experimental tensile test data
    data_path = Path('data/experimental/surfalex_tensile_test_data.csv')
    data = read_tensile_test_data(data_path)

    # Use the stress/strain data for the first orientation (zero degrees)
    ori_idx = 0

    exp_tensile_test = TensileTest(
        eng_stress=data[ori_idx]['eng_stress'] * 1e6,
        eng_strain=data[ori_idx]['eng_strain'],
    )

    with Path('./data/lm_opt/damask_material_params.yml').open() as handle:
        damask_params = YAML().load(handle)

    # Which parameters should be optimised, and where are they in the parameters dict:
    fitting_params = [
        FittingParameter('h0_slipslip', [400e6], ['phase', 0, 'keys'], 0.1),
        FittingParameter('tausat_slip', [95e6], ['phase', 0, 'keys'], 0.1),
    ]

    # Tell the fitter what module and function to use to write the simulation input files:
    inputs_writer = {
        'module': 'damask_parse',
        'function': 'write_material_config',
        'function_args': {
            'material': '<<PARAMETERS>>',
            'dir_path': '<<SIM_DIR>>',
            'part_paths': {
                'Texture': './texture.config',
                'Microstructure': './microstructure.config',
            },
        },
    }

    # Generate a new LMFitter object:
    lm_fitter = LMFitter(
        exp_tensile_test,
        damask_params,
        fitting_params,
        inputs_writer,
        base_sim_dir=Path('./data/lm_opt/base_damask_sim'),
        sim_dir=Path('./data/lm_opt/sims/'),
        initial_damping=[2, 1, 0.5],
    )

    lm_fitter.generate_simulation_inputs()
    lm_fitter.to_json_file(LM_FIT_JSON_NAME)


def iterate():
    'Continue the optimisation process'

    print('{}.iterate'.format(__name__))

    lm_fitter = LMFitter.from_json_file(LM_FIT_JSON_NAME)

    # Load results of previous set of simulations:
    com = (lm_fitter.opt_index - 1) * lm_fitter.sims_per_iteration
    sim_dir_range = range(com, com + lm_fitter.sims_per_iteration)
    tensile_tests = []
    for i in sim_dir_range:
        table_path = lm_fitter.sim_dir.joinpath(i, 'postProc', 'geom_load.txt')
        table_data = read_table(table_path)
        tensile_tests.append(
            TensileTest(
                true_stress=table_data['Mises(Cauchy)'],
                true_strain=table_data['Mises(ln(V))']
            )
        )

    lm_fitter.add_simulated_tensile_tests(tensile_tests)
    lm_fitter.generate_simulation_inputs()
    lm_fitter.to_json_file(LM_FIT_JSON_NAME)


if __name__ == '__main__':
    if sys.argv[1] == 'set-up':
        set_up()
    elif sys.argv[1] == 'iterate':
        iterate()
