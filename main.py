# -*- coding: utf-8 -*-

"""
STALICLA
PBPM Module (2)
This script calls the PBPM module and integrates its functions with the user arguments.
Author: Gabriela Martinez - airamgabriela17@gmail.com
"""

from src.PBPM import *


def main():

    """ Main function to execute the PBPM protocol """

    # ------ Arguments parser ------ #

    args = create_arg_parser()

    # ------ Directories setting ------ #

    directory = os.path.dirname(os.path.abspath(__file__))

    default_folder_output = 'analysis'
    int_folder_output = 'data/processed'
    folder_logs = 'analysis'

    dir_default_folder_output = os.path.join(directory, default_folder_output)
    dir_int_folder_output = os.path.join(directory, int_folder_output)
    dir_logs = os.path.join(directory, folder_logs)

    if args['path'] is not None:
        dir_pbpm = args['path']
    else:
        dir_pbpm = dir_default_folder_output

    # ------ Filenames and parameters ------ #

    csv_name = os.path.join(dir_int_folder_output, 'base-matrix-{}.csv'.format(strftime("%Y-%m-%d_%H꞉%m꞉%S", gmtime())))
    b_csv_name = os.path.join(dir_pbpm, 'binary-matrix-{}.csv'.format(strftime("%Y-%m-%d_%H꞉%m꞉%S", gmtime())))
    n_csv_name = os.path.join(dir_pbpm, 'numerical-matrix-{}.csv'.format(strftime("%Y-%m-%d_%H꞉%m꞉%S", gmtime())))
    nn_csv_name = os.path.join(dir_pbpm, 'normalized-matrix-{}.csv'.format(strftime("%Y-%m-%d_%H꞉%m꞉%S", gmtime())))

    command = "# Script generated with the following parameters --- " + \
              " inputFile: " + args['inputFile'] + \
              " pathwaysDirectory: " + args['pathwaysDirectory'] + \
              "--pathways: " + none_to_str(args['pathway']) + \
              "--genes: " + none_to_str(args['gene']) + \
              "--patients: " + none_to_str(args['patient']) + \
              "--consequences: " + none_to_str(args['csq']) + \
              "--pli: " + none_to_str(args['pli_gt']) + \
              "--pr: " + none_to_str(args['pr_g']) + \
              "--af: " + none_to_str(args['af_lt']) + \
              "--polyphen: " + none_to_str(args['pph2']) + \
              "--mpc: " + none_to_str(args['mpc_gt']) + \
              "--adjusted_consequences: " + none_to_str(args['adj_csq']) + \
              "--matrix: " + none_to_str(args['matrix']) + \
              "--intermediate_matrix: " + bool_to_str(args['im']) + \
              "--append: " + none_to_str(args['append']) + \
              "--path: " + none_to_str(args['path'])

    # ------ Data processing and filtering ------ #

    # Input files reading, filtering and cleaning
    patients = drop_na(read_mutations(args['inputFile']))
    pathways = drop_na(read_pathways(args['pathwaysDirectory'], args['pathway']))
    genes_to_norm = genes_to_normalize_pathways(pathways)

    # Joining data: patient mutations and pathways
    joint_data = join_input_data(patients, pathways)

    # Restrict analysis to specific patients
    data = filter_patients(joint_data, args['patient'])

    # Restrict analysis to specific genes
    data = filter_genes(data, args['gene'])

    # Restrict analysis to specific consequences
    data = filter_consequences(data, args['csq'])

    # Restrict analysis to specific pLI
    data = filter_pli(data, args['pli_gt'])

    # Restrict analysis to specific pRecessive
    data = filter_precessive(data, args['pr_g'])

    # Restrict analysis to specific max_control_AF
    data = filter_af(data, args['af_lt'])

    # Restrict analysis to specific values of pph2 predictions
    data = filter_polyphen(data, args['pph2'])

    # Restrict analysis to specific values of MPC
    data = filter_mpc(data, args['mpc_gt'])

    # Restrict analysis to specific values of adjusted consequence
    data = filter_adj_consequence(data, args['adj_csq'])

    # ------ Generate intermediate matrix ------ #

    # Create intermediate matrix
    int_matrix = create_base_matrix(data)

    # Download intermediate matrix if desired
    download_matrix(args['im'], int_matrix, command, csv_name)

    if args['append'] is None:

        # ------ Generate PBPM matrix ------ #

        # Compute desired PBPM matrix type and download a copy of it
        pbpm = create_pbpm(args['matrix'], data, genes_to_norm)
        download_pbpm(args['matrix'], pbpm, command, b_csv_name, n_csv_name, nn_csv_name)

    else:

        ap_b_csv_name = os.path.join(dir_pbpm, 'appended-binary-matrix-{}.csv'.format(strftime("%Y-%m-%d_%H꞉%m꞉%S", gmtime())))
        ap_n_csv_name = os.path.join(dir_pbpm, 'appended-numerical-matrix-{}.csv'.format(strftime("%Y-%m-%d_%H꞉%m꞉%S", gmtime())))
        ap_nn_csv_name = os.path.join(dir_pbpm, 'appended-normalized-matrix-{}.csv'.format(strftime("%Y-%m-%d_%H꞉%m꞉%S", gmtime())))

        # ------ Generate intermediate appended matrix ------ #

        # Create an appended matrix if desired
        int_matrix_appended = append_base_matrices(int_matrix, command, args['append'])

        # Return appended matrix (for internal processing)
        int_matrix_appended_to_return = return_file(int_matrix_appended)

        # Process appended matrix for end user and export a copy of it
        int_matrix_appended_path = edit_lines(args['append'], int_matrix_appended_to_return, dir_int_folder_output)

        # ------ Generate new PBPM appended matrix ------ #

        # Recompute PBPM and download it
        appended_pbpm = append_to_pbpm(int_matrix_appended_path, pathways, args['matrix'], genes_to_norm)
        download_pbpm(args['matrix'], appended_pbpm, command, ap_b_csv_name, ap_n_csv_name, ap_nn_csv_name)

    # ------ Logging of execution ------ #

    logging_info(args, dir_logs, dir_pbpm, dir_int_folder_output)


if __name__ == "__main__":
    main()
