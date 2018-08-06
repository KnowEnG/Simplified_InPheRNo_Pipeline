"""
lanier4@illinois.edu

"""
import os

def InPheRNo_simple(run_parameters):
    """ wrapper for calling InPheRNo from the command line """

    valid_keys = ['id', 'it', 'ie', 'igp', 'od', 'onp', 'ons']

    k_plus = ' -'
    s_plus = ' '

    # relative path from the "test" directory             ---     Fix Me
    src_dir = os.path.abspath('../src')
    script_filename = os.path.join(src_dir, 'InPheRNo_simplified.py')
    run_command = 'python3 ' +  script_filename

    # fails here unless run_parameters contain all valid keys
    for k in valid_keys:
        run_command += k_plus + k + s_plus + run_parameters[k]

    os.system(run_command)


SELECT = {"InPheRNo_simple": InPheRNo_simple}


def main():
    """
    This is the main function to perform sample clustering
    """
    import sys
    from knpackage.toolbox import get_run_directory_and_file
    from knpackage.toolbox import get_run_parameters

    run_directory, run_file = get_run_directory_and_file(sys.argv)
    run_parameters = get_run_parameters(run_directory, run_file)

    SELECT[run_parameters["method"]](run_parameters)


if __name__ == "__main__":
    main()
