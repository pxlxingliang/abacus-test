import argparse, importlib, os, sys
from abacustest.lib_model.model_args import MODEL_ARGS


def _import_model_class(model_info):
    """Import model class from MODEL_ARGS info, with error handling."""
    try:
        module = importlib.import_module(f"abacustest.lib_model.{model_info['file']}")
        return getattr(module, model_info['class_name'])
    except ImportError as e:
        print(f"Error: Failed to import model '{model_info['file']}'")
        print(f"Import error: {e}")
        print("This may be due to missing dependencies. Please install the required packages.")
        sys.exit(1)
    except AttributeError as e:
        print(f"Error: Failed to find class '{model_info['class_name']}' in model '{model_info['file']}'")
        print(f"Attribute error: {e}")
        sys.exit(1)


def ModelArgs(parser, list_all_models=False, model_subcommand=None):  
    subparser = parser.add_subparsers(dest="model")
    
    if list_all_models:
        model_names = MODEL_ARGS.keys()
        
        for model_name in model_names:
            model_info = MODEL_ARGS[model_name]
            iparser = subparser.add_parser(model_name, help=model_info["description"])
            iparser.description = model_info["description"]
    else:
        if model_subcommand is None:
            return
        
        if model_subcommand not in MODEL_ARGS:
            print(f"Error: unknown model '{model_subcommand}'")
            print(f"Available models: {', '.join(sorted(MODEL_ARGS.keys()))}")
            sys.exit(1)
        
        model_info = MODEL_ARGS[model_subcommand]
        model_class = _import_model_class(model_info)
        
        iparser = subparser.add_parser(model_subcommand, help=model_info["description"])
        iparser.description = model_info["description"]
        
        if model_class.HAS_PREPARE_POST_COMMAND:
            subsubparser = iparser.add_subparsers(dest="modelcommand")
            iparser_prepare = subsubparser.add_parser("prepare", help="Prepare the model")
            iparser_postprocess = subsubparser.add_parser("post", help = "Post-process the model")
            model_class.prepare_args(iparser_prepare)
            model_class.postprocess_args(iparser_postprocess)

        model_class.add_args(iparser)
    
    parser.description = "Prepare and post-process the specified model"

    return parser

def RunModel(param):
    model_info = MODEL_ARGS[param.model]
    model_class = _import_model_class(model_info)
    imodel = model_class()
    
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