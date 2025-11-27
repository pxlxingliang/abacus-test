import argparse,importlib,glob,os
from abacustest.lib_model.model import Model

def ReadAllModels():
    return {imodel.model_name():imodel for imodel in Model.__subclasses__()} 

def ModelArgs(parser):  
    allmodels = [os.path.basename(i)[:-3] for i in glob.glob(os.path.join(os.path.dirname(__file__),"lib_model/model_*.py"))]
    for module in allmodels:
        importlib.import_module(f"abacustest.lib_model.{module}")   
    subparser = parser.add_subparsers(dest="model")
    
    models = sorted(Model.__subclasses__(), key=lambda x: x.model_name().lower())
    
    for imodel in models:
        model_name = imodel.model_name()
        iparser = subparser.add_parser(model_name, help=imodel.description())
        iparser.description = imodel.description()
        
        if imodel.HAS_PREPARE_POST_COMMAND:
            subsubparser = iparser.add_subparsers(dest="modelcommand")
            iparser_prepare = subsubparser.add_parser("prepare", help="Prepare the model")
            iparser_postprocess = subsubparser.add_parser("post", help = "Post-process the model")
            imodel.prepare_args(iparser_prepare)
            imodel.postprocess_args(iparser_postprocess)

        imodel.add_args(iparser)
    parser.description = "Prepare and post-process the specified model"

    return parser

def RunModel(param):
    allmodels = ReadAllModels()
    imodel = allmodels[param.model]()
    
    print(f"Model: {param.model}")
    if not imodel.HAS_PREPARE_POST_COMMAND:
        imodel.run(param)
    else:
        if param.modelcommand == "prepare":
            imodel.run_prepare(param)
        elif param.modelcommand == "post":
            imodel.run_postprocess(param)
        else:
            imodel.run(param)

def main():
    parser = argparse.ArgumentParser()
    param = ModelArgs(parser).parse_args()
    RunModel(param)
    
if __name__ == "__main__":
    main()