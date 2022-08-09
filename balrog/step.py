# Making a first bit of code for Eastlake 

import eastlake

from eastlake.step import Step
from eastlake.des_files import read_pizza_cutter_yaml

########################## Line Width Limit ###################################

class BalrogGalSimRunner(Step):
    """
    Pipeline step which runs galsim
    The config attribute is a little different here, since it is updated when
    running GalSim
    """

    def __init__(
        self, config, base_dir, name="Balrog", logger=None, verbosity=0, 
        log_file=None):
        
        super().__init__(
            config, base_dir, name=name, logger=logger, verbosity=verbosity,
            log_file=log_file)


        
    def execute(self, stash, new_params=None, except_abort=False, verbosity=1.,
                log_file=None, comm=None):

        if comm is not None:
            rank = comm.Get_rank()
        else:
            rank = 0

        if new_params is not None:
            galsim.config.UpdateConfig(self.config, new_params)
           
        if rank == 0:
            self.logger.debug(
                "Process config dict: \n%s", pprint.pformat(config))
          
        if self.name not in stash:
            stash[self.name] = {}

        # Get the tilename
        # TODO: Need to fix for Balrog, need to get the tilename in the way
        # Balrog does - megan and Brian fix this 
        stash["tilenames"] = [config["output"]["tilename"]]

        # TODO: Brian and megan need to run Balrog here 
        #galsim.config.Process(config, self.logger, except_abort=except_abort)

        self.update_stash(config, stash)

        # Return status and stash
        return 0, stash
 


    def update_stash(self, config, stash):

        # Update the stash with information on image files etc. required by
        # following steps

        # Get the output type and number of files
        # TODO: Brian and megan need to have Balrog say which bands and tiles
        #bands = config["output"]["bands"]
        #nbands = len(bands)
        #tilenames = stash["tilenames"]
        #tilename = tilenames[0]
        assert len(tilenames) == 1

        self.logger.error(
            "Simulated tile %s in bands %s" % (tilename, str(bands)))
        stash["nbands"] = nbands # set a list of stings fo bands TODO
        stash["bands"] = bands # See above TODO

        # TODO Brian and megan need to set this 
        stash["desrun"] = desrun # The string like y3v02 but for Y6

        # Need to set this - is the input data 
        imsim_data = os.environ['IMSIM_DATA'] 
        stash["imsim_data"] = imsim_data

        #base_dir = self.base_dir # keep if want

        #put for bands in band: TODO

        # TODO: What is this line?
        stash.set_input_pizza_cutter_yaml(
                        read_pizza_cutter_yaml(imsim_data, desrun, tilename, band),
                        tilename,
                        band)

        # truth
        # This is where we write out the info for what Balrog did
        with stash.update_output_pizza_cutter_yaml(tilename, band) as pyml:
            for i in range(len(pyml["src_info"])):
                fname = pyml["src_info"][i]["image_path"]
                if fname.endswith(".fz"):
                    fname = fname[:-3]
                    
        # TODO: MEGAN WHAT LEVEL IS THIS CODE AT??? 
        # WHERE DID THIS CODE COME FROM? 
        # Src_info is a dictionary of SE images i is the SE images
        # Writes the images uncompressed here 	
        pyml["src_info"][i]["image_path"] = fname
        pyml["src_info"][i]["image_ext"] = 0

        pyml["src_info"][i]["bmask_path"] = fname
        pyml["src_info"][i]["bmask_ext"] = 1

        pyml["src_info"][i]["weight_path"] = fname
        pyml["src_info"][i]["weight_ext"] = 2

        # Delete these for now
        #truth_files = [
          #  get_truth_from_image_file(src["image_path"], tilename)
           # for src in pyml["src_info"]
        #]
        stash.set_filepaths("truth_files", truth_files, tilename, band=band)

        # also get tile center
        tile_center = get_tile_center(
            stash.get_input_pizza_cutter_yaml(tilename, bands[0])["image_path"])
        stash.set_tile_info_quantity("tile_center", tile_center, tilename)
    
    
        return 



    @classmethod
    def from_config_file(cls, config_file, logger=None):
        all_config = galsim.config.ReadConfig(config_file, None, logger)
        assert len(all_config) == 1
        return cls(all_config[0], logger=logger)



    def set_base_dir(self, base_dir):

        # This is where the output gets written - tell Balrog to work here 
        self.base_dir = base_dir

        # Update the output directory. - dont do this - MB
        #self.config['output']['dir'] = base_dir




eastlake.register_pipeline_step("balrog", BalrogGalSimRunner, is_galsim=False)




















