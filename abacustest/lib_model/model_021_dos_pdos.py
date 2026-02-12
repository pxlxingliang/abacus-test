from ..model import Model
import os, traceback
from abacustest.lib_model.comm_dos import DOSData, PDOSData


class DOSPDOSModel(Model):
    @staticmethod
    def model_name(): # type: ignore
        '''
        Name of the model, which will be used as the subcommand
        '''
        return "dos-pdos"
    
    @staticmethod
    def description(): # type: ignore
        '''
        Description of the model
        '''
        return "Postprocess the DOS and PDOS data from ABACUS calculations"
    
    @staticmethod
    def add_args(parser):
        '''
        Add arguments for the model
        '''
        parser.description = "Postprocess the DOS and PDOS calculation. Will analyze and plot DOS and PDOS data."
        parser.add_argument('-j', '--job', default=None, nargs="+", help='the path of abacus inputs (required)')
        parser.add_argument("--range", default=[-10,10], type=float, help="The energy range for plots, default is -10,10 eV", nargs=2)
        parser.add_argument("--atom-index", default=None, type=int, help="Atom index to analyze (1-based). If not specified, all atoms will be analyzed.")
        parser.add_argument("--plot-type", default="species", type=str, choices=["species", "shell", "orbital", "atom"], help="Plot type: species, shell, orbital, or atom")
        parser.add_argument("--suffix", default=None, type=str, help="Suffix for output files. If set, files will be named DOS_{suffix}.dat/png and PDOS_{suffix}.dat/png")
        parser.add_argument("--no-save-data", action="store_true", help="Do not save data files (default is to save)")
        parser.add_argument("--no-save-plot", action="store_true", help="Do not save plot files (default is to save)")
        return parser
    
    def run(self, params):
        '''
        Run the postprocess process
        '''
        if params.job is None:
            print("Error: Please specify job directory with -j/--job")
            return 1
        
        jobs = params.job
        if len(jobs) == 0:
            print("No valid job found")
            return 1
            
        PostDOSPDOS(jobs, 
                   energy_range=params.range,
                   atom_index=params.atom_index,
                   plot_type=params.plot_type,
                   suffix=params.suffix,
                   save_data=not params.no_save_data,
                   save_plot=not params.no_save_plot).run()
        print("DOS and PDOS analysis completed.")
        return 0


class PostDOSPDOS:
    def __init__(self, jobs, energy_range=None, 
                 atom_index=None, plot_type="species", suffix=None, save_data=False, save_plot=False):
        self.jobs = jobs
        self.energy_range = [-10, 10] if energy_range is None else energy_range
        self.atom_index = atom_index
        self.plot_type = plot_type
        self.suffix = suffix
        self.save_data = save_data
        self.save_plot = save_plot
                    
    def run(self):
        for job in self.jobs:
            print(f"Processing job: {job}")
            try:
                self.process_single_job(job)
            except Exception as e:
                print(f"Error processing job {job}: {e}")
                traceback.print_exc()
    
    def _get_filename(self, job, base_name, ext):
        """Generate filename with optional suffix."""
        if self.suffix:
            return os.path.join(job, f"{base_name}_{self.suffix}.{ext}")
        else:
            return os.path.join(job, f"{base_name}.{ext}")
    
    def process_single_job(self, job):
        # Read INPUT file from job directory
        inputfile = os.path.join(job, "INPUT")
        if not os.path.isfile(inputfile):
            print(f"Warning: Cannot find INPUT file in {job}")
            return
            
        # Read DOS and PDOS data
        try:
            dos_data = DOSData.ReadFromAbacusJob(job)
            pdos_data = PDOSData.ReadFromAbacusJob(job)
        except Exception as e:
            print(f"Warning: Could not read DOS/PDOS data from {job}: {e}")
            return
        
        # Plot DOS and PDOS according to parameters
        if self.save_plot:
            self.plot_dos(dos_data, job)
            
            if self.plot_type == "species":
                self.plot_species_pdos(pdos_data, job)
            elif self.plot_type == "shell":
                self.plot_shell_pdos(pdos_data, job)
            elif self.plot_type == "orbital":
                self.plot_orbital_pdos(pdos_data, job)
            elif self.plot_type == "atom":
                self.plot_atom_pdos(pdos_data, job)
            
        # Save data if requested
        if self.save_data:
            self.save_dos_data(dos_data, job)
            self.save_pdos_data(pdos_data, job)
    
    def plot_dos(self, dos_data, job):
        """Plot DOS"""
        try:
            plot_filename = self._get_filename(job, "DOS", "png")
            dos_data.plot_dos(
                emin=self.energy_range[0],
                emax=self.energy_range[1],
                title="Density of States",
                fname=plot_filename
            )
            print(f"Saved DOS plot to {plot_filename}")
        except Exception as e:
            print(f"Error plotting DOS: {e}")
    
    def plot_species_pdos(self, pdos_data, job):
        """Plot PDOS for all species"""
        try:
            species = pdos_data.get_species()
            if len(species) == 0:
                return
            
            plot_filename = self._get_filename(job, "PDOS", "png")
            
            pdos_data.plot_species_pdos(
                emin=self.energy_range[0], 
                emax=self.energy_range[1], 
                pdos_fig_name=plot_filename
            )
            print(f"Saved species PDOS plot to {plot_filename}")
        except Exception as e:
            print(f"Error plotting species PDOS: {e}")
    
    def plot_shell_pdos(self, pdos_data, job):
        """Plot PDOS for shells of all species"""
        try:
            species = pdos_data.get_species()
            if len(species) == 0:
                return
            
            plot_filename = self._get_filename(job, "PDOS", "png")
            
            pdos_data.plot_species_shell_pdos(
                emin=self.energy_range[0], 
                emax=self.energy_range[1], 
                pdos_fig_name=plot_filename
            )
            print(f"Saved shell PDOS plot to {plot_filename}")
        except Exception as e:
            print(f"Error plotting shell PDOS: {e}")
    
    def plot_orbital_pdos(self, pdos_data, job):
        """Plot PDOS for orbitals of all species"""
        try:
            species = pdos_data.get_species()
            if len(species) == 0:
                return
            
            plot_filename = self._get_filename(job, "PDOS", "png")
            
            pdos_data.plot_species_orbital_pdos(
                emin=self.energy_range[0], 
                emax=self.energy_range[1], 
                pdos_fig_name=plot_filename
            )
            print(f"Saved orbital PDOS plot to {plot_filename}")
        except Exception as e:
            print(f"Error plotting orbital PDOS: {e}")
    
    def plot_atom_pdos(self, pdos_data, job):
        """Plot PDOS for specific atoms"""
        try:
            # Get unique atom indices from PDOS data (1-based indexing from ABACUS)
            unique_atom_indices = sorted(set([orb['atom_index'] for orb in pdos_data.projected_dos]))
            
            if self.atom_index is not None:
                # Use user-provided 1-based index directly
                if self.atom_index not in unique_atom_indices:
                    print(f"Warning: Atom index {self.atom_index} not found in PDOS data")
                    return
                atom_indices = [self.atom_index]
            else:
                # Use first 3 atoms (already 1-based)
                atom_indices = unique_atom_indices[:3]
            
            if len(atom_indices) == 0:
                return
            
            plot_filename = self._get_filename(job, "PDOS", "png")
            
            pdos_data.plot_atoms_pdos(
                atom_indices=atom_indices,
                emin=self.energy_range[0], 
                emax=self.energy_range[1], 
                pdos_fig_name=plot_filename
            )
            print(f"Saved atom PDOS plot to {plot_filename}")
        except Exception as e:
            print(f"Error plotting atom PDOS: {e}")
            traceback.print_exc()
    
    def save_pdos_data(self, pdos_data, job):
        """Save PDOS data to files"""
        try:
            # Get unique atom indices from PDOS data (1-based indexing from ABACUS)
            unique_atom_indices = sorted(set([orb['atom_index'] for orb in pdos_data.projected_dos]))
            pdos_file = self._get_filename(job, "PDOS", "dat")
            
            if self.plot_type == "species":
                pdos_data.write_species_pdos(pdos_dat_file=pdos_file)
                print(f"Saved species PDOS data to {pdos_file}")
            elif self.plot_type == "shell":
                pdos_data.write_species_shell_pdos(pdos_dat_file=pdos_file)
                print(f"Saved shell PDOS data to {pdos_file}")
            elif self.plot_type == "orbital":
                pdos_data.write_species_orbital_pdos(pdos_dat_file=pdos_file)
                print(f"Saved orbital PDOS data to {pdos_file}")
            elif self.plot_type == "atom":
                if self.atom_index is not None:
                    if self.atom_index not in unique_atom_indices:
                        print(f"Warning: Atom index {self.atom_index} not found in PDOS data")
                        return
                    pdos_data.write_atoms_pdos(atom_indices=[self.atom_index], pdos_dat_file=pdos_file)
                    print(f"Saved atom PDOS data to {pdos_file}")
                else:
                    atom_indices = unique_atom_indices[:3]
                    if len(atom_indices) > 0:
                        pdos_data.write_atoms_pdos(atom_indices=atom_indices, pdos_dat_file=pdos_file)
                        print(f"Saved atom PDOS data to {pdos_file}")
        except Exception as e:
            print(f"Error saving PDOS data: {e}")
            traceback.print_exc()

    def save_dos_data(self, dos_data, job):
        """Save DOS data to file"""
        try:
            dos_filename = self._get_filename(job, "DOS", "dat")
            dos_data.write_dos(dos_filename)
            print(f"Saved DOS data to {dos_filename}")
        except Exception as e:
            print(f"Error saving DOS data: {e}")
            traceback.print_exc()
