# File: read_namelist.py
import f90nml


#warning: this function is not used anywhere
def readNamelist(pathNamelist):
    """
    Reads a Fortran namelist file and returns its contents as a dictionary.

    Parameters:
        namelist_path (str): Path to the namelist file.

    Returns:
        dict: Dictionary containing the namelist parameters.
    """
    try:
        nml = f90nml.read(pathNamelist)
        return nml
    
    except Exception as e:
        print(f"Error reading namelist: {e}")
        return None