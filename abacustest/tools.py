import argparse, importlib, sys
from abacustest.lib_tools.tools_args import TOOLS_ARGS


def _import_tool_class(tool_info):
    """Import tool class from TOOLS_ARGS info, with error handling."""
    try:
        module = importlib.import_module(f"abacustest.lib_tools.{tool_info['file']}")
        return getattr(module, tool_info['class_name'])
    except ImportError as e:
        print(f"Error: Failed to import tool '{tool_info['file']}'")
        print(f"Import error: {e}")
        print("This may be due to missing dependencies. Please install the required packages.")
        sys.exit(1)
    except AttributeError as e:
        print(f"Error: Failed to find class '{tool_info['class_name']}' in tool '{tool_info['file']}'")
        print(f"Attribute error: {e}")
        sys.exit(1)


def ToolsArgs(parser, list_all_tools=False, tool_subcommand=None):
    subparser = parser.add_subparsers(dest="tool")

    if list_all_tools:
        for tool_name in TOOLS_ARGS:
            tool_info = TOOLS_ARGS[tool_name]
            iparser = subparser.add_parser(tool_name, help=tool_info["description"])
            iparser.description = tool_info["description"]
    else:
        if tool_subcommand is None:
            return

        if tool_subcommand not in TOOLS_ARGS:
            print(f"Error: unknown tool '{tool_subcommand}'")
            print(f"Available tools: {', '.join(sorted(TOOLS_ARGS.keys()))}")
            sys.exit(1)

        tool_info = TOOLS_ARGS[tool_subcommand]
        tool_class = _import_tool_class(tool_info)

        iparser = subparser.add_parser(tool_subcommand, help=tool_info["description"])
        iparser.description = tool_info["description"]
        tool_class.add_args(iparser)

    parser.description = "Miscellaneous utility tools"

    return parser


def RunTool(param):
    tool_info = TOOLS_ARGS[param.tool]
    tool_class = _import_tool_class(tool_info)
    itool = tool_class()

    print(f"Tool: {param.tool}")
    itool.run(param)


def main():
    parser = argparse.ArgumentParser()
    param = ToolsArgs(parser).parse_args()
    RunTool(param)


if __name__ == "__main__":
    main()
