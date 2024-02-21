import numpy as np

default_config = {
    "GLOBAL": {"PROJECT": "DPGEN"},
    "FORCE_EVAL": {
        "METHOD": "QS",
        "STRESS_TENSOR": "ANALYTICAL",
        "DFT": {
            "BASIS_SET_FILE_NAME": "./cp2k_basis_pp_file/BASIS_MOLOPT",
            "POTENTIAL_FILE_NAME": "./cp2k_basis_pp_file/GTH_POTENTIALS",
            "CHARGE": 0,
            "UKS": "F",
            "MULTIPLICITY": 1,
            "MGRID": {"CUTOFF": 400, "REL_CUTOFF": 50, "NGRIDS": 4},
            "QS": {"EPS_DEFAULT": "1.0E-12"},
            "SCF": {"SCF_GUESS": "ATOMIC", "EPS_SCF": "1.0E-6", "MAX_SCF": 50},
            "XC": {"XC_FUNCTIONAL": {"_": "PBE"}},
        },
        "SUBSYS": {
            "CELL": {"A": "10 .0 .0", "B": ".0 10 .0", "C": ".0 .0 10"},
            "COORD": {"@include": "coord.xyz"},
            "KIND": [
                {
                    "_": "O",
                    "POTENTIAL": "GTH-PBE-q6",
                    "BASIS_SET": "DZVP-MOLOPT-SR-GTH"
                }, {
                    "_": "H",
                    "POTENTIAL": "GTH-PBE-q1",
                    "BASIS_SET": "DZVP-MOLOPT-SR-GTH"
                }, {
                    "_": "Pt",
                    "POTENTIAL": "GTH-PBE-q10",
                    "BASIS_SET": "DZVP-A5-Q10-323-MOL-T1-DERIVED_SET-1"
                }, {
                    "_": "Ag",
                    "POTENTIAL": "GTH-PBE-q11",
                    "BASIS_SET": "DZVP-MOLOPT-SR-GTH"
                }, {
                    "_": "Na",
                    "POTENTIAL": "GTH-PBE-q9",
                    "BASIS_SET": "DZVP-MOLOPT-SR-GTH"
                }, {
                    "_": "K",
                    "POTENTIAL": "GTH-PBE-q9",
                    "BASIS_SET": "DZVP-MOLOPT-SR-GTH"
                }, {
                    "_": "Li",
                    "POTENTIAL": "GTH-PBE-q3",
                    "BASIS_SET": "DZVP-MOLOPT-SR-GTH"
                }, {
                    "_": "C",
                    "POTENTIAL": "GTH-PBE-q4",
                    "BASIS_SET": "DZVP-MOLOPT-SR-GTH"
                }, {
                    "_": "N",
                    "POTENTIAL": "GTH-PBE-q5",
                    "BASIS_SET": "DZVP-MOLOPT-SR-GTH"
                }, {
                    "_": "Cl",
                    "POTENTIAL": "GTH-PBE-q7",
                    "BASIS_SET": "DZVP-MOLOPT-SR-GTH"
                }, {
                    "_": "F",
                    "POTENTIAL": "GTH-PBE-q7",
                    "BASIS_SET": "DZVP-MOLOPT-SR-GTH"
                }, {
                    "_": "Mg",
                    "POTENTIAL": "GTH-PBE-q10",
                    "BASIS_SET": "DZVP-MOLOPT-SR-GTH"
                }, {
                    "_": "Al",
                    "POTENTIAL": "GTH-PBE-q3",
                    "BASIS_SET": "DZVP-MOLOPT-SR-GTH"
                }
            ]
        },
        "PRINT": {"FORCES": {"_": "ON"}, "STRESS_TENSOR": {"_": "ON"}},
    },
}


def update_dict(old_d, update_d):
    """A method to recursive update dict
    :old_d: old dictionary
    :update_d: some update value written in dictionary form.
    """
    import collections.abc

    for k, v in update_d.items():
        if (
            k in old_d
            and isinstance(old_d[k], dict)
            and isinstance(update_d[k], collections.abc.Mapping)
        ):
            update_dict(old_d[k], update_d[k])
        else:
            old_d[k] = update_d[k]


def iterdict(input_dict, out_list, loop_idx):
    """ 
    recursively generate a list of strings for further 
    print out CP2K input file

    Args:
        input_dict: dictionary for CP2K input parameters
        out_list: list of strings for printing
        loop_idx: record of loop levels in recursion
    Return:
        out_list
    """
    if len(out_list) == 0:
        out_list.append("\n")
    start_idx = len(out_list) - loop_idx - 2
    for k, v in input_dict.items():
        k = str(k)  # cast key into string
        #if value is dictionary
        if isinstance(v, dict):
            out_list.insert(-1 - loop_idx, "&" + k)
            out_list.insert(-1 - loop_idx, "&END " + k)
            iterdict(v, out_list, loop_idx + 1)
        # if value is list
        elif isinstance(v, list):
            if isinstance(v[0], dict):
                for _v in v:
                    out_list.insert(-1 - loop_idx, "&" + k)
                    out_list.insert(-1 - loop_idx, "&END " + k)
                    iterdict(_v, out_list, loop_idx + 1)
            else:
                for _v in v:
                    _v = str(_v)
                    out_list.insert(-1 - loop_idx, k + " " + _v)
        # if value is other type, e.g., int/float/str
        else:
            v = str(v)
            if k == "_":
                out_list[start_idx] = out_list[start_idx] + " " + v
            else:
                out_list.insert(-1 - loop_idx, k + " " + v)
    return out_list


def make_cp2k_input(sys_data, fp_params):
    # covert cell to cell string
    cell = sys_data["cells"][0]
    cell = np.reshape(cell, [3, 3])
    cell_a = np.array2string(cell[0, :])
    cell_a = cell_a[1:-1]
    cell_b = np.array2string(cell[1, :])
    cell_b = cell_b[1:-1]
    cell_c = np.array2string(cell[2, :])
    cell_c = cell_c[1:-1]

    # get update from user
    user_config = fp_params
    # get update from cell
    cell_config = {
        "FORCE_EVAL": {"SUBSYS": {"CELL": {"A": cell_a, "B": cell_b, "C": cell_c}}}
    }
    update_dict(default_config, user_config)
    update_dict(default_config, cell_config)
    # output list
    input_str = iterdict(default_config, out_list=["\n"], loop_idx=0)
    string = "\n".join(input_str)
    string = string.strip("\n")
    return string


def make_cp2k_xyz(sys_data):
    # get structral information
    atom_names = sys_data["atom_names"]
    atom_types = sys_data["atom_types"]

    # write coordinate to xyz file used by cp2k input
    coord_list = sys_data["coords"][0]
    u = np.array(atom_names)
    atom_list = u[atom_types]
    x = "\n"
    for kind, coord in zip(atom_list, coord_list):
        x += str(kind) + " " + str(coord[:])[1:-1] + "\n"
    return x


def make_cp2k_input_from_external(sys_data, exinput_path):
    # read the input content as string
    with open(exinput_path) as f:
        exinput = f.readlines()

    # find the ABC cell string
    for line_idx, line in enumerate(exinput):
        if "ABC" in line:
            delete_cell_idx = line_idx
            delete_cell_line = line

    # remove the useless CELL line
    exinput.remove(delete_cell_line)

    # insert the cell information
    # covert cell to cell string
    cell = sys_data["cells"][0]
    cell = np.reshape(cell, [3, 3])
    cell_a = np.array2string(cell[0, :])
    cell_a = cell_a[1:-1]
    cell_b = np.array2string(cell[1, :])
    cell_b = cell_b[1:-1]
    cell_c = np.array2string(cell[2, :])
    cell_c = cell_c[1:-1]

    exinput.insert(delete_cell_idx, "A  " + cell_a + "\n")
    exinput.insert(delete_cell_idx + 1, "B  " + cell_b + "\n")
    exinput.insert(delete_cell_idx + 2, "C  " + cell_c + "\n")

    return "".join(exinput)
