#!/usr/bin/env python3
"""
Script to generate test cases for abacus collectdata.

Usage:
    python generate_test.py <path_to_abacus_example> [-o output_file.py]

Rules for test generation:
1. INPUT: Only verify it's a dict, and check calculation/esolver_type/kspacing values
2. Lists with >3 elements: Only verify length and first/middle/last elements
3. band_plot: Skip (path-dependent)
"""

import sys
import os
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from abacustest.lib_collectdata.collectdata import RESULT

INPUT_SPECIAL_KEYS = ["calculation", "esolver_type", "kspacing"]
SKIP_METRICS = ["band_plot"]


def is_path_dependent(val):
    """Check if a value is a path that depends on cwd"""
    if not isinstance(val, str):
        return False
    path_indicators = ['/', '\\', '.png', '.txt', '.dat', '.xlsx', '.cube', '.html']
    return any(ind in val for ind in path_indicators) or Path(val).is_absolute()


def is_long_list(val):
    """Check if value is a list with more than 3 elements OR is deeply nested"""
    if not isinstance(val, list):
        return False
    if len(val) > 3:
        return True
    if len(val) > 0 and isinstance(val[0], list):
        depth = get_list_depth(val)
        if depth >= 2:
            return True
    return False


def generate_assertion(name, val):
    """Generate assertion line for a metric"""
    if val is None:
        return f"        # {name}: None"

    if name == "INPUT":
        return generate_input_assertion(val)
    elif name in SKIP_METRICS:
        return f"        # {name}: skipped"
    elif is_long_list(val):
        return generate_long_list_assertion(name, val)
    elif isinstance(val, bool):
        return f'        self.assertEqual(ab["{name}"], {val})'
    elif isinstance(val, int):
        return f'        self.assertEqual(ab["{name}"], {val})'
    elif isinstance(val, float):
        return f'        self.assertAlmostEqual(ab["{name}"], {val}, places=3)'
    elif isinstance(val, str):
        if is_path_dependent(val):
            return f'        # {name}: "{val}" (path-dependent, skipped)'
        return f'        self.assertEqual(ab["{name}"], "{val}")'
    elif isinstance(val, list) and len(val) == 0:
        return f'        self.assertEqual(ab["{name}"], [])'
    elif isinstance(val, list):
        if is_long_list(val):
            return generate_long_list_assertion(name, val)
        else:
            # Check if this is a flat list of floats
            all_float = all(isinstance(x, float) for x in val)
            if all_float:
                lines = []
                lines.append(f"        self.assertEqual(len(ab[\"{name}\"]), {len(val)})")
                for i, v in enumerate(val):
                    lines.append(f'        self.assertAlmostEqual(ab["{name}"][{i}], {v}, places=3)')
                return "\n".join(lines)
            # Check if this is a 2D list of floats
            if all(isinstance(row, list) for row in val) and all(isinstance(x, float) for row in val for x in row):
                lines = []
                lines.append(f"        self.assertEqual(len(ab[\"{name}\"]), {len(val)})")
                for i, row in enumerate(val):
                    lines.append(f"        self.assertEqual(len(ab[\"{name}\"][{i}]), {len(row)})")
                    for j, v in enumerate(row):
                        lines.append(f'        self.assertAlmostEqual(ab["{name}"][{i}][{j}], {v}, places=3)')
                return "\n".join(lines)
            first = val[0]
            if isinstance(first, (int, float)):
                return f'        self.assertEqual(ab["{name}"], {val})'
            elif isinstance(first, str):
                return f'        self.assertEqual(ab["{name}"], {val})'
            else:
                return f'        self.assertEqual(ab["{name}"], {val})'
    elif isinstance(val, dict):
        return f'        self.assertEqual(ab["{name}"], {val})'
    else:
        return f'        # {name}: {type(val).__name__} = {val}'


def generate_input_assertion(val):
    """Generate assertion for INPUT metric"""
    lines = []
    lines.append("        input_dict = ab[\"INPUT\"]")
    lines.append("        self.assertIsInstance(input_dict, dict)")
    for key in INPUT_SPECIAL_KEYS:
        if key in val:
            v = val[key]
            if isinstance(v, str):
                lines.append(f'        self.assertEqual(input_dict["{key}"], "{v}")')
            else:
                lines.append(f'        self.assertEqual(input_dict["{key}"], {v})')
    return "\n".join(lines)


def get_list_depth(val):
    """Get nesting depth of a list"""
    if not isinstance(val, list):
        return 0
    if len(val) == 0:
        return 1
    return 1 + get_list_depth(val[0])


def flatten_list(val):
    """Flatten a nested list into 1D"""
    if not isinstance(val, list):
        return [val]
    result = []
    for item in val:
        result.extend(flatten_list(item))
    return result


def generate_long_list_assertion(name, val):
    """Generate assertion for long lists (>3 elements)"""
    length = len(val)
    lines = []

    first = val[0]
    depth = get_list_depth(val)

    if depth == 1:
        length = len(val)
        mid_index = length // 2
        middle = val[mid_index]
        lines.append(f"        self.assertEqual(len(ab[\"{name}\"]), {length})")
        if isinstance(first, (int, float)):
            lines.append(f"        self.assertAlmostEqual(ab[\"{name}\"][0], {first}, places=3)")
            lines.append(f"        self.assertAlmostEqual(ab[\"{name}\"][len(ab[\"{name}\"]) // 2], {middle}, places=3)")
            lines.append(f"        self.assertAlmostEqual(ab[\"{name}\"][-1], {val[-1]}, places=3)")
        elif isinstance(first, str):
            lines.append(f'        self.assertEqual(ab[\"{name}\"][0], "{first}")')
            lines.append(f'        self.assertEqual(ab[\"{name}\"][len(ab[\"{name}\"]) // 2], "{middle}")')
            lines.append(f'        self.assertEqual(ab[\"{name}\"][-1], "{val[-1]}")')
        else:
            lines.append(f"        # ab[\"{name}\"][0] = {first}")
    elif depth == 2:
        lines.append(f"        self.assertEqual(len(ab[\"{name}\"]), {length})")
        lines.append(f"        self.assertEqual(len(ab[\"{name}\"][0]), {len(first)})")
        flat = [item for sublist in val for item in sublist]
        lines.append(f"        {name}_flat = [item for sublist in ab[\"{name}\"] for item in sublist]")
        lines.append(f"        self.assertEqual(len({name}_flat), {len(flat)})")
        lines.append(f"        self.assertAlmostEqual({name}_flat[0], {flat[0]}, places=3)")
        lines.append(f"        self.assertAlmostEqual({name}_flat[len({name}_flat) // 2], {flat[len(flat) // 2]}, places=3)")
        lines.append(f"        self.assertAlmostEqual({name}_flat[-1], {flat[-1]}, places=3)")
    elif depth >= 3:
        lines.append(f"        self.assertEqual(len(ab[\"{name}\"]), {length})")
        lines.append(f"        self.assertEqual(len(ab[\"{name}\"][0]), {len(first)})")
        lines.append(f"        self.assertEqual(len(ab[\"{name}\"][0][0]), {len(first[0])})")
        flat = [item for sublist in val for subsublist in sublist for item in subsublist]
        lines.append(f"        {name}_flat = [item for sublist in ab[\"{name}\"] for subsublist in sublist for item in subsublist]")
        lines.append(f"        self.assertEqual(len({name}_flat), {len(flat)})")
        lines.append(f"        self.assertAlmostEqual({name}_flat[0], {flat[0]}, places=3)")
        lines.append(f"        self.assertAlmostEqual({name}_flat[len({name}_flat) // 2], {flat[len(flat) // 2]}, places=3)")
        lines.append(f"        self.assertAlmostEqual({name}_flat[-1], {flat[-1]}, places=3)")
    else:
        lines.append(f"        # ab[\"{name}\"][0] = {first}")

    return "\n".join(lines)


def generate_test(path, output_file=None):
    """Generate test code from abacus example"""
    ab = RESULT(path=path, fmt="abacus")

    from abacustest.lib_collectdata.abacus.abacus import Abacus
    all_names = sorted(Abacus.AllMethod().keys())

    not_none = []
    is_none = []

    for name in all_names:
        try:
            val = ab[name]
            if val is not None:
                not_none.append((name, val))
            else:
                is_none.append(name)
        except Exception:
            is_none.append(name)

    lines = []
    lines.append("import unittest")
    lines.append("from pathlib import Path")
    lines.append("")
    lines.append("from abacustest.lib_collectdata.collectdata import RESULT")
    lines.append("")
    lines.append("")
    lines.append("class TestAbacusCollectdata(unittest.TestCase):")
    lines.append('    """Test Abacus result collection from SCF calculation"""')
    lines.append("")
    lines.append("    def test_abacus_scf(self):")
    lines.append('        """Test all extractable metrics from abacus-scf example"""')
    lines.append(f'        ab = RESULT(path=Path(__file__).parent / "{path}", fmt="abacus")')
    lines.append("")

    for name, val in not_none:
        lines.append(generate_assertion(name, val))

    lines.append("")
    lines.append("        # Metrics that returned None for this example:")
    for name in is_none:
        lines.append(f"        # {name}: None")

    lines.append("")
    lines.append("")
    lines.append("if __name__ == '__main__':")
    lines.append("    unittest.main()")

    result = "\n".join(lines)

    if output_file:
        with open(output_file, 'w') as f:
            f.write(result)
        print(f"Test file written to: {output_file}")
    else:
        print(result)

    return not_none, is_none


def main():
    if len(sys.argv) < 2:
        print("Usage: python generate_test.py <path_to_abacus_example> [-o output_file.py]")
        print("Example: python generate_test.py ./abacus-scf")
        print("Example: python generate_test.py ./abacus-scf -o test_output.py")
        sys.exit(1)

    path = sys.argv[1]
    output_file = None

    if len(sys.argv) >= 4 and sys.argv[2] == "-o":
        output_file = sys.argv[3]

    generate_test(path, output_file)


if __name__ == "__main__":
    main()
