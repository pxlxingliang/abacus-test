"""
ABACUS Numerical Atomic Orbitals (NAO) file handling.
"""

from typing import Dict, List, Any, Optional

import numpy as np


class AbacusNAO:
    """
    Class for handling ABACUS Numerical Atomic Orbitals (NAO) files.

    Attributes:
        element: Element symbol
        energy_cutoff: Energy cutoff in Ry
        radius: Radius cutoff in atomic units (Bohr)
        lmax: Maximum angular momentum quantum number
        l_orbs: List of number of orbitals for each angular momentum
                l_orbs[0] = number of s orbitals (l=0)
                l_orbs[1] = number of p orbitals (l=1)
                etc.
        mesh: Number of grid points
        dr: Grid spacing
        orbs: List of orbital data, each containing:
              - l: angular momentum quantum number
              - n: principal quantum number
              - data: radial wavefunction values
    """

    def __init__(
        self,
        element: str,
        energy_cutoff: float,
        radius: float,
        lmax: int,
        l_orbs: List[int],
        mesh: int,
        dr: float,
        orbs: Optional[List[Dict[str, Any]]] = None,
    ) -> None:
        self.element = element
        self.energy_cutoff = energy_cutoff
        self.radius = radius
        self.lmax = lmax
        self.l_orbs = l_orbs
        self.mesh = mesh
        self.dr = dr
        self.orbs = orbs if orbs is not None else []

    @staticmethod
    def read_from_file(nao_file: str) -> "AbacusNAO":
        """
        Read a NAO file and return an AbacusNAO instance.

        Args:
            nao_file: Path to the NAO file

        Returns:
            AbacusNAO instance
        """
        with open(nao_file, "r") as f:
            lines = f.readlines()

        summary_end_line = 0
        for linenum, line in enumerate(lines):
            if "SUMMARY  END" in line:
                summary_end_line = linenum
                break

        if summary_end_line == 0:
            raise ValueError(f"Could not find 'SUMMARY  END' in {nao_file}")

        summary = lines[1:summary_end_line]
        element = summary[0].split()[-1]
        energy_cutoff = float(summary[1].split()[-1])
        radius = float(summary[2].split()[-1])
        lmax = int(summary[3].split()[-1])

        l_orbs = []
        for line in summary[4:]:
            line_stripped = line.strip()
            if line_stripped.startswith("Number of"):
                l_orbs.append(int(line.split()[-1]))

        mesh = None
        dr = None
        for i in range(summary_end_line + 1, min(summary_end_line + 5, len(lines))):
            line = lines[i]
            if "Mesh" in line:
                mesh = int(line.split()[-1])
            elif "dr" in line:
                dr = float(line.split()[-1])

        if mesh is None or dr is None:
            raise ValueError(f"Could not find Mesh or dr in {nao_file}")

        orb_l_n_lines = []
        for i in range(summary_end_line + 1, len(lines)):
            line = lines[i]
            if "Type" in line and "L" in line and "N" in line:
                orb_l_n_lines.append(i + 1)

        orbs = []
        for orb_idx, orb_l_n_line_idx in enumerate(orb_l_n_lines):
            tokens = lines[orb_l_n_line_idx].strip().split()
            l = int(tokens[1])
            n = int(tokens[2])

            data_start_line = orb_l_n_line_idx + 1
            if orb_idx < len(orb_l_n_lines) - 1:
                data_end_line = orb_l_n_lines[orb_idx + 1] - 1
            else:
                data_end_line = len(lines)

            orb_data = []
            for line in lines[data_start_line:data_end_line]:
                if not line.strip():
                    continue
                try:
                    for num in line.split():
                        orb_data.append(float(num))
                except ValueError:
                    continue

            orbs.append({"l": l, "n": n, "data": np.array(orb_data)})

        return AbacusNAO(
            element=element,
            energy_cutoff=energy_cutoff,
            radius=radius,
            lmax=lmax,
            l_orbs=l_orbs,
            mesh=mesh,
            dr=dr,
            orbs=orbs,
        )

    @property
    def basis_num(self) -> int:
        """
        Total number of basis functions.
        Assume the number of orbitals for s, p, d, f, g, h are 1, 3, 5, 7, 9, 11 respectively.
        """
        nbas = 0
        for i, norb in enumerate(self.l_orbs):
            nbas += (2 * i + 1) * norb

        return nbas

    def get_orbital_label(self, l: int, n: int) -> str:
        """
        Get the orbital label (e.g., '1s', '2p', '3d', '4f').

        Args:
            l: Angular momentum quantum number
            n: Principal quantum number

        Returns:
            Orbital label string
        """
        labels = {0: "s", 1: "p", 2: "d", 3: "f", 4: "g", 5: "h"}
        return f"{n + 1}{labels.get(l, '?')}"

    def __repr__(self) -> str:
        return (
            f"AbacusNAO(element={self.element}, energy_cutoff={self.energy_cutoff}, "
            f"radius={self.radius}, lmax={self.lmax}, l_orbs={self.l_orbs}, "
            f"mesh={self.mesh}, dr={self.dr})"
        )
