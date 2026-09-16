"""
bsse: Prepare the STRU files needed for BSSE (Basis Set Superposition Error) correction.

Given an input STRU file and the atom indices of one fragment (the other fragment is
the complement), this tool generates 5 directories:

    BSSE_AB          : the original full structure (fragment A + fragment B)
    BSSE_A_ghostB    : all atoms, fragment B atoms turned into empty atoms (H_empty, O_empty, ...)
    BSSE_B_ghostA    : all atoms, fragment A atoms turned into empty atoms
    BSSE_A           : fragment A atoms only
    BSSE_B           : fragment B atoms only

Each directory contains the corresponding 'STRU' file. If an ABACUS INPUT file is
given or an 'INPUT' file exists in the current directory, it is copied into every
BSSE directory as well (otherwise no INPUT is copied).

Empty (ghost) atoms keep the mass, pseudopotential and numerical orbital of their
original element, so they only contribute basis functions without any ionic potential.

The BSSE energy can then be evaluated as:
    BSSE = E(A_ghostB) + E(B_ghostA) - E(A) - E(B)

NOTE:
- Atom indices given by the user are 1-based (the first atom in the STRU file is 1).
- Only ABACUS STRU files are supported as input.
"""
from abacustest.lib_tools.tool import Tool
from abacustest.lib_prepare.stru import AbacusSTRU
import copy
import os
import shutil


def _parse_indices(index_str):
    """Parse a 1-based atom index string such as '1,3-5' into a sorted list of integers.

    Args:
        index_str (str): Comma separated list of integers or integer ranges (e.g. '1,2,3-5').

    Returns:
        list[int]: Sorted list of 1-based atom indices.

    Raises:
        ValueError: If the string contains invalid syntax or non-positive indices.
    """
    indices = []
    for part in index_str.replace(" ", "").split(","):
        if not part:
            continue
        if "-" in part:
            start, end = part.split("-", 1)
            try:
                start, end = int(start), int(end)
            except ValueError:
                raise ValueError(f"Invalid index range: '{part}'")
            if start <= 0 or end <= 0 or start > end:
                raise ValueError(f"Invalid index range: '{part}' (indices are 1-based)")
            indices.extend(range(start, end + 1))
        else:
            try:
                idx = int(part)
            except ValueError:
                raise ValueError(f"Invalid index: '{part}'")
            if idx <= 0:
                raise ValueError(f"Invalid index: '{part}' (indices are 1-based)")
            indices.append(idx)
    return sorted(set(indices))


def _make_ghost(stru, ghost_indices):
    """Return a copy of the structure where the given atoms are turned into empty atoms.

    Empty atoms keep all the properties (element, mass, pp, orb, mag, ...) of the
    original atom but their label is changed to '<element>_empty'. The returned
    structure is sorted so that atoms of the same type are grouped, which is required
    by the ABACUS STRU format.

    Args:
        stru (AbacusSTRU): The original structure.
        ghost_indices (list[int]): 0-based atom indices to be turned into empty atoms.

    Returns:
        AbacusSTRU: A new structure with the ghost atoms relabeled.
    """
    atoms = []
    for i, atom in enumerate(stru.atoms):
        new_atom = copy.deepcopy(atom)
        if i in ghost_indices:
            if new_atom.element is None:
                raise ValueError(
                    f"Cannot infer element for atom index {i + 1} with label "
                    f"'{new_atom.label}', so its empty label cannot be determined."
                )
            new_atom.label = f"{new_atom.element}_empty"
        atoms.append(new_atom)

    new_stru = AbacusSTRU(
        cell=stru.cell,
        atoms=atoms,
        dpks=stru.dpks,
        metadata=copy.deepcopy(stru.metadata),
    )
    new_stru.sort()
    return new_stru


class BsseTool(Tool):
    @staticmethod
    def add_args(parser):
        parser.description = (
            "Prepare the STRU files needed for BSSE (counterpoise) correction. "
            "Given a STRU file and the atom indices of one fragment, generate 5 BSSE_* "
            "directories each containing a STRU file: the full structure, the full "
            "structure with each fragment turned into empty atoms, and the two isolated "
            "fragments. An ABACUS INPUT file is copied into each directory when given "
            "or when present in the current directory."
        )
        parser.add_argument(
            "-i", "--input", default="STRU",
            help="Input ABACUS STRU file."
        )
        parser.add_argument(
            "-o", "--output-dir", dest="output_dir", default=None,
            help="Directory in which the BSSE_* directories are created. "
                 "Default: the current directory."
        )
        parser.add_argument(
            "--input-file", dest="input_file", default=None,
            help="ABACUS INPUT file to copy into each BSSE directory. "
                 "Default: the 'INPUT' file in the current directory if present, "
                 "otherwise no INPUT is copied."
        )
        parser.add_argument(
            "--frag1", default=None,
            help="1-based atom indices of fragment A, comma separated and/or with ranges, "
                 "e.g. '1,2,3-5'. If only one fragment is given, the remaining atoms form "
                 "the other fragment."
        )
        parser.add_argument(
            "--frag2", default=None,
            help="1-based atom indices of fragment B. Optional; if both --frag1 and --frag2 "
                 "are given, they must be disjoint and cover all atoms."
        )
        return parser

    def run(self, params):
        input_file = params.input
        if not os.path.isfile(input_file):
            print(f"Error: input file '{input_file}' does not exist.")
            return

        stru = AbacusSTRU.read(input_file, fmt="stru")
        if stru is None:
            print(f"Error: failed to read STRU file '{input_file}'.")
            return

        natoms = stru.natoms
        if natoms == 0:
            print("Error: the input structure has no atoms.")
            return

        try:
            frag1 = _parse_indices(params.frag1) if params.frag1 else None
            frag2 = _parse_indices(params.frag2) if params.frag2 else None
        except ValueError as e:
            print(f"Error: {e}")
            return

        if frag1 is None and frag2 is None:
            print("Error: at least one of --frag1 or --frag2 must be provided.")
            return

        if frag1 is None:
            frag1 = sorted(set(range(1, natoms + 1)) - set(frag2))
        if frag2 is None:
            frag2 = sorted(set(range(1, natoms + 1)) - set(frag1))

        for name, idxs in (("frag1", frag1), ("frag2", frag2)):
            for i in idxs:
                if not (1 <= i <= natoms):
                    print(f"Error: {name} contains invalid index {i}; "
                          f"valid 1-based indices are 1-{natoms}.")
                    return

        if not frag1:
            print("Error: fragment A (frag1) has no atoms.")
            return
        if not frag2:
            print("Error: fragment B (frag2) has no atoms.")
            return

        overlap = sorted(set(frag1) & set(frag2))
        if overlap:
            print(f"Error: overlapping indices between frag1 and frag2: {overlap}")
            return
        if len(frag1) + len(frag2) != natoms:
            print(f"Error: frag1 ({len(frag1)}) + frag2 ({len(frag2)}) must cover "
                  f"all {natoms} atoms.")
            return

        frag1_0 = [i - 1 for i in frag1]
        frag2_0 = [i - 1 for i in frag2]

        output_dir = params.output_dir if params.output_dir else os.getcwd()
        os.makedirs(output_dir, exist_ok=True)

        input_file_src = params.input_file
        if input_file_src is None:
            cwd_input = os.path.join(os.getcwd(), "INPUT")
            if os.path.isfile(cwd_input):
                input_file_src = cwd_input
        if input_file_src is not None and not os.path.isfile(input_file_src):
            print(f"Error: INPUT file '{input_file_src}' does not exist.")
            return

        print(f"Input structure: {input_file} ({natoms} atoms)")
        print("Atom indices (1-based) and their fragment assignment:")
        for i, atom in enumerate(stru.atoms):
            frag = "A" if i in frag1_0 else "B"
            print(f"  {i + 1:>4}: {atom.label:<12} -> fragment {frag}")

        stru_ghost_b = _make_ghost(stru, frag2_0)
        stru_ghost_a = _make_ghost(stru, frag1_0)
        stru_a = stru.create_subset(frag1_0)
        stru_b = stru.create_subset(frag2_0)
        stru_a.sort()
        stru_b.sort()

        outputs = [
            ("BSSE_AB", stru, "E(A+B), full complex", "copy"),
            ("BSSE_A_ghostB", stru_ghost_b,
             "E(A in AB basis), fragment B as empty/ghost", None),
            ("BSSE_B_ghostA", stru_ghost_a,
             "E(B in AB basis), fragment A as empty/ghost", None),
            ("BSSE_A", stru_a, "E(A), fragment A alone", None),
            ("BSSE_B", stru_b, "E(B), fragment B alone", None),
        ]

        for name, s, desc, mode in outputs:
            out_dir = os.path.join(output_dir, name)
            os.makedirs(out_dir, exist_ok=True)
            stru_path = os.path.join(out_dir, "STRU")
            if mode == "copy":
                shutil.copy2(input_file, stru_path)
            else:
                if not s.write(stru_path, fmt="stru"):
                    print(f"Error: failed to write '{stru_path}'.")
                    return
            if input_file_src is not None:
                shutil.copy2(input_file_src, os.path.join(out_dir, "INPUT"))
            print(f"Written {stru_path} ({s.natoms} atoms)  # {desc}")

        if input_file_src is not None:
            print(f"Copied INPUT file into each BSSE directory.")
        else:
            print("No INPUT file found; no INPUT copied into the BSSE directories.")

        print("\nBSSE = E(A_ghostB) + E(B_ghostA) - E(A) - E(B)")
        print("Remember to adjust the 'ntype' parameter in the INPUT file for each "
              "calculation.")
