class termcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def print_info(str):
    print termcolors.OKBLUE + str + termcolors.ENDC

def print_warning(str):
    print termcolors.WARNING + str + termcolors.ENDC

def print_error(str):
    print termcolors.FAIL + str + termcolors.ENDC
