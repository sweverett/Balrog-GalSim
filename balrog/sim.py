def make_images(conf, prep):
    """
    parameters
    -----------
    conf: dict
        The configuration dictionary.

    prep: obj
        A Preparator object

        It contains the info about the  null weight files,
        psfs, etc.

        It can be used like a dict, e.g. it has prep['nwgint_flist'] available

    """
    maker=ImageMaker(conf, prep)
    maker.go()

class ImageMaker(dict):
    """
    class to run galsim and make the new images

    parameters
    ----------
    config: dict
        The configuration dictionary
    """
    def __init__(self, conf, prep):
        self.update(conf)
        self.prep=prep

    def go(self):
        """
        make the images
        """
        pass
