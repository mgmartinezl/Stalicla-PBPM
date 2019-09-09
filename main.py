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

    dir_default_folder_output = os.path.join(directory, default_folder_output)
    dir_int_folder_output = os.path.join(directory, int_folder_output)

    if args['path'] is not None:
        dir_pbpm = args['path']
        dir_logs = args['path']
    else:
        dir_pbpm = dir_default_folder_output
        dir_logs = dir_default_folder_output

    # ------ Filenames and parameters ------ #

    csv_name = os.path.join(dir_int_folder_output, 'base-matrix-{}.csv'.format(local_time.strftime("%Y-%m-%d_%H꞉%m꞉%S")))

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
    patients = format_df(drop_na(read_mutations(args['inputFile'])))
    pathways = format_df(drop_na(read_pathways(args['pathwaysDirectory'], args['pathway'])))
    genes_to_norm = genes_to_normalize_pathways(pathways)

    # Joining data: patient mutations and pathways
    joint_data = join_input_data(patients, pathways)

    # Restrict analysis to specific patients
    data = filter_patients(joint_data, args['patient'])

    # ------ Provisional PBPM ------ #

    # Compute a first PBPM only filtered for specific patients or pathways
    pbpm_1 = create_pbpm(args['matrix'], data, genes_to_norm)

    # Make all columns null for the first PBPM
    pbpm_1 = make_cols_null(pbpm_1)

    # ------ Provisional PBPM filtering ------ #

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
    data = filter_PolyPhen(data, args['pph2'])

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

        # ------ Generate final PBPM ------ #

        # Compute a second PBPM only for filtering parameters
        pbpm_2 = create_pbpm(args['matrix'], data, genes_to_norm)

        # Combine matrices and download
        pbpm = combine_matrices(pbpm_1, pbpm_2)
        print(pbpm_name(args['matrix']))
        download_pbpm(pbpm, command, os.path.join(dir_pbpm, pbpm_name(args['matrix'])))

    else:

        # ------ Generate intermediate appended matrix ------ #

        # Create an appended matrix if desired
        int_matrix_appended = append_base_matrices(int_matrix, command, args['append'])

        # Return appended matrix (for internal processing)
        int_matrix_appended_to_return = return_file(int_matrix_appended)

        # Process appended matrix for end user and export a copy of it
        int_matrix_appended_path = edit_and_export_file(args['append'], int_matrix_appended_to_return,
                                                        dir_int_folder_output)

        # ------ Generate final appended PBPM ------ #

        # Compute a second PBPM only for filtering parameters
        appended_joint_data = join_input_data(format_appended_file(int_matrix_appended_path), pathways)
        appended_pbpm = create_pbpm(args['matrix'], appended_joint_data, genes_to_norm)

        # Combine matrices and download
        appended_pbpm = combine_matrices(pbpm_1, appended_pbpm)
        download_pbpm(appended_pbpm, command, os.path.join(dir_pbpm, pbpm_name(args['matrix'])))

    # ------ Logging of execution ------ #

    logging_info(args, os.path.join(dir_logs, log_name(args['matrix'])), dir_pbpm, dir_int_folder_output)


if __name__ == "__main__":
    main()
