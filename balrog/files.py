import yaml

def read_yaml(fname):
    """
    read a yaml file
    """
    with open(fname) as fobj:
        data=yaml.load(fobj)

    return data
