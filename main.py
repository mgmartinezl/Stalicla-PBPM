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
    # Arguments parser
    args = create_arg_parser()

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

    csv_name = os.path.join(dir_int_folder_output, 'base-matrix-{}.csv'.format(strftime("%Y-%m-%d_%H꞉%m꞉%S", gmtime())))
    b_csv_name = os.path.join(dir_pbpm, 'binary-matrix-{}.csv'.format(strftime("%Y-%m-%d_%H꞉%m꞉%S", gmtime())))
    n_csv_name = os.path.join(dir_pbpm, 'numerical-matrix-{}.csv'.format(strftime("%Y-%m-%d_%H꞉%m꞉%S", gmtime())))
    nn_csv_name = os.path.join(dir_pbpm, 'normalized-matrix-{}.csv'.format(strftime("%Y-%m-%d_%H꞉%m꞉%S", gmtime())))

    command = "# Script generated with the following parameters --- " + \
              " inputFile: " + args['inputFile'] + \
              " pathwaysDirectory: " + args['pathwaysDirectory'] + \
              "-pathway: " + none_to_str(args['pathway']) + \
              "-gene: " + none_to_str(args['gene']) + \
              "-patient: " + none_to_str(args['patient']) + \
              "-mutation: " + none_to_str(args['csq']) + \
              "-pli (greater/equal than): " + none_to_str(args['pli_gt']) + \
              "-af (less/equal than): " + none_to_str(args['af_lt']) + \
              "-pph2: " + none_to_str(args['pph2']) + \
              "-mpc (greater/equal than): " + none_to_str(args['mpc_gt']) + \
              "-adj_csq: " + none_to_str(args['adj_csq']) + \
              "-matrix: " + none_to_str(args['matrix']) + \
              "-append: " + none_to_str(args['append']) + \
              "-path: " + none_to_str(args['path'])

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

    # Compute desired PBPM matrix type and download a copy of it
    pbpm = create_pbpm(args['matrix'], data, genes_to_norm)
    download_pbpm(args['matrix'], pbpm, command, b_csv_name, n_csv_name, nn_csv_name)

    # Create raw base matrix
    base_matrix = create_base_matrix(data)

    # Download base matrix if desired
    base_matrix_downloaded = download_matrix(args['im'], base_matrix, command, csv_name)

    # Create an appended matrix if desired
    base_matrix_appended = append_base_matrices(base_matrix_downloaded, command, args['append'])

    # Return appended matrix (for internal processing)
    base_matrix_appended_to_return = return_file(base_matrix_appended)

    # Process appended matrix for end user
    edit_lines(args['append'], base_matrix_appended_to_return, dir_int_folder_output)

    logging_info(args, dir_logs, dir_pbpm, dir_int_folder_output)


if __name__ == "__main__":
    main()
