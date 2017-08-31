def make_images(prep):
    """
    parameters
    -----------
    prep: a preparator objects
        prep is a preparator object.  It contains the config info and can be
        used as a dict, in addition to methods

        e.g. it has prep['nwgint_flist'] available


    """
    maker=ImageMaker(prep)
    maker.go()

class ImageMaker(dict):
    """
    class to run galsim and make the new images
    """
    def __init__(self, prep):
        self.update(prep)

    def go(self):
        """
        make the images
        """
