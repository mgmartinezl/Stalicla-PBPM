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

    directory = os.path.dirname(os.path.abspath(__file__))
    folder_output = 'reports'
    folder_logs = 'analysis'
    dir_output = os.path.join(directory, folder_output)
    dir_logs = os.path.join(directory, folder_logs)

    csv_name = os.path.join(folder_output, 'base-matrix-{}.csv'.format(strftime("%Y-%m-%d_%H꞉%m꞉%S", gmtime())))
    b_csv_name = os.path.join(folder_output, 'binary-matrix-{}.csv'.format(strftime("%Y-%m-%d_%H꞉%m꞉%S", gmtime())))
    n_csv_name = os.path.join(folder_output, 'numerical-matrix-{}.csv'.format(strftime("%Y-%m-%d_%H꞉%m꞉%S", gmtime())))
    nn_csv_name = os.path.join(folder_output, 'normalized-matrix-{}.csv'.format(strftime("%Y-%m-%d_%H꞉%m꞉%S", gmtime())))

    # Arguments parser
    args = create_arg_parser()

    command = "# Script generated with the following parameters --- " + \
              " inputFile: " + args['inputFile'] + \
              " pathwaysDirectory: " + args['pathwaysDirectory'] + \
              "-pathway: " + none_to_str(args['pathway']) + \
              "-gene: " + none_to_str(args['gene']) + \
              "-patient: " + none_to_str(args['patient']) + \
              "-mutation: " + none_to_str(args['mutation']) + \
              "-pli (greater/equal than): " + none_to_str(args['pli_gt']) + \
              "-af (less than): " + none_to_str(args['af_lt']) + \
              "-pph2: " + none_to_str(args['pph2']) + \
              "-mpc (greater/equal than): " + none_to_str(args['mpc_gt']) + \
              "-adj (consequence): " + none_to_str(args['adj']) + \
              "-matrix: " + none_to_str(args['matrix']) + \
              "-r " + none_to_str(args['r']) + \
              "-append " + none_to_str(args['append'])

    # Input files reading, filtering and cleaning
    patients = drop_na(read_mutations(args['inputFile']))
    pathways = drop_na(read_pathways(args['pathwaysDirectory'], args['pathway']))
    genes_to_norm = genes_to_normalize_pathways(pathways)

    # Joining data: patient mutations and pathways
    joint_data = join_input_data(patients, pathways)

    # Restrict analysis to specific patients
    data = analyze_patients(joint_data, args['patient'])

    # Restrict analysis to specific genes
    data = analyze_genes(data, args['gene'])

    # Restrict analysis to specific mutations
    data = analyze_mutations(data, args['mutation'])

    # Restrict analysis to specific pLI
    data = analyze_pli(data, args['pli_gt'])

    # Restrict analysis to specific max_control_AF
    data = analyze_max_control_af(data, args['af_lt'])

    # Restrict analysis to specific values of pph2 predictions
    data = analyze_pph2_prediction(data, args['pph2'])

    # Restrict analysis to specific values of MPC
    data = analyze_mpc(data, args['mpc_gt'])

    # Restrict analysis to specific values of adjusted consequence
    data = analyze_adj_consequence(data, args['adj'])

    # Compute desired PBPM matrix type and download a copy of it
    pbpm = create_pbpm(args['matrix'], data, genes_to_norm)

    download_pbpm(args['matrix'], pbpm, command, b_csv_name, n_csv_name, nn_csv_name)

    # Create raw base matrix
    base_matrix = create_base_matrix(data)

    # Download base matrix if desired
    base_matrix_downloaded = download_matrix(args['r'], base_matrix, command, csv_name)

    # Create an appended matrix if desired
    base_matrix_appended = append_base_matrices(base_matrix_downloaded, command, args['append'])

    # Return appended matrix (for internal processing)
    base_matrix_appended_to_return = return_file(base_matrix_appended)

    # Process appended matrix for end user
    edit_lines(args['append'], base_matrix_appended_to_return, dir_output)

    logging_info(args, dir_logs, dir_output)


if __name__ == "__main__":
    main()
