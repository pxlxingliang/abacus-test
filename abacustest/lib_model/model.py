class Model:
    '''
    Base class for all models
    
    Some specific workflow will be implemented in the subclass.
    Such as EOS, BAND, etc.
    
    The model should have the following methods:
    '''
    HAS_PREPARE_POST_COMMAND = True


    @staticmethod
    def model_name():
        '''
        Name of the model, which will be used as the subcommand
        '''
        return "Model"
    
    @staticmethod
    def description():
        '''
        Description of the model
        '''
        return ""
    
    @staticmethod
    def add_args(parser):
        '''
        Add arguments for the model
        The arguments can not be command, model, modelcommand 
        '''
        pass
    
    def run(self,params):
        '''
        Parse the parameters and run the model
        '''
        pass
    
    @staticmethod
    def prepare_args(parser):
        '''
        Add arguments for the prepare subcommand
        The arguments can not be command, model, modelcommand '''
        pass
    
    def run_prepare(self,params):
        '''
        Parse the parameters and run the prepare process.
        Usually, this step will generate the input files for abacustest submit.
        '''
        pass
    
    @staticmethod
    def postprocess_args(parser):
        '''
        Add arguments for the postprocess subcommand
        The arguments can not be command, model, modelcommand'''
        pass

    def run_postprocess(self,params):
        '''
        Parse the parameters and run the postprocess process'''
        pass