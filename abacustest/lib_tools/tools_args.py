"""
TOOLS_ARGS - Tools registry for abacustest

This dictionary defines all available tools in abacustest.
Each tool is identified by a unique key (tool_name) that is used as the subcommand.

Format of each entry:
    "tool_name": {
        "description": "Brief description of the tool (shown in help)",
        "file": "Filename of the tool module (without .py extension)",
        "class_name": "Name of the Tool subclass in the module"
    }

Adding a new tool:
1. Create a new file in abacustest/lib_tools/ (e.g., tool_xxx_name.py)
2. Define a class that inherits from Tool (from abacustest.lib_tools.tool import Tool)
3. Implement the required methods:
   - add_args(parser): Add command-line arguments
   - run(params): Implement tool functionality
4. Add an entry to TOOLS_ARGS with:
   - Key: The tool name (used as subcommand)
   - description: Brief description
   - file: Module filename (without .py)
   - class_name: The class name (must match the class definition)

Note:
- The tool_name in the key must be unique and should be lowercase with hyphens if needed
- The file should follow the naming pattern: tool_###_Name.py (e.g., tool_001_supercell.py)
- The class_name must exactly match the class defined in the module
"""

TOOLS_ARGS = {
    "supercell": {
        "description": "Extend the unit cell to supercell",
        "file": "tool_001_supercell",
        "class_name": "SuperCellTool",
    },
    "vasp2abacus": {
        "description": "Transform VASP input files to ABACUS input files.",
        "file": "tool_002_vasp2abacus",
        "class_name": "Vasp2AbacusTool",
    },
    "conv": {
        "description": "Convert structure files between formats (POSCAR, CIF, ABACUS STRU)",
        "file": "tool_003_conv",
        "class_name": "StructureConvertTool",
    },
}
