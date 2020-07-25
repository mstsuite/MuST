"""ASE interface for MST module implemented in MuST package available at
 https://github.com/mstsuite/MuST"""

import numpy as np
from ase.calculators.calculator import FileIOCalculator, SCFError
import os
import subprocess
import glob
from ase.units import Bohr, Rydberg
from ase_must.default_params import defaults
from ase.data import atomic_numbers

magmoms = {'Fe': 2.1,
           'Co': 1.4,
           'Ni': 0.6}


def write_positions_input(atoms, method):
    """ Function that writes the positions data input file based on
    atoms object and selected calculation method"""
    with open('position.dat', 'w') as filehandle:
        filehandle.write(str(1.0) + '\n\n')

        for i in range(3):
            filehandle.write('%s\n' % str(atoms.get_cell()[i] / Bohr)[1:-1])
        filehandle.write('\n')

        if method == 3:

            for site in atoms.info['CPA']:
                sitestring = 'CPA  %s' % str(atoms[site['index']].position
                                             / Bohr)[1:-1]

                for key in site.keys():
                    if key == 'index':
                        pass
                    else:
                        sitestring += '  %s %s' % (key, str(site[key]))
                sitestring += '\n'
                filehandle.write(sitestring)

        else:
            for index in range(len(atoms)):
                filehandle.write('%s %s\n'
                                 % (atoms[index].symbol,
                                    str(atoms[index].position / Bohr)[1:-1]))


def write_atomic_pot_input(symbol, nspins, moment, xc, niter, mp):
    """
    Function to write input file for generating
    atomic potential using 'newa' command
    Parameters
    ----------
    symbol: str
        Chemical symbol of the element.
    nspins: int
        Number of spins.
    moment: float
        Magnetic moment.
    xc: int
        ex-cor type (1=vb-hedin,2=vosko).
    niter: int
        Maximum number of SCF iterations.
    mp: float
        SCF mixing parameter
    """
    title = symbol + ' Atomic Potential'
    output_file = symbol + '_a_out'
    pot_file = symbol + '_a_pot'
    z = atomic_numbers[symbol]

    if moment == 0. and nspins == 2 and symbol in ['Fe', 'Co', 'Ni']:
        moment = magmoms[symbol]

    space = 20 * ' '
    contents = [('Title', title),
                (output_file,
                 'Output file name. If blank, data will show on screen'),
                (z, 'Atomic number'),
                (moment, 'Magnetic moment'),
                (nspins, 'Number of spins'),
                (xc, 'Exchange-correlation type (1=vb-hedin,2=vosko)'),
                (niter, 'Number of Iterations'),
                (mp, 'Mixing parameter'),
                (pot_file, 'Output potential file')]

    with open(symbol + '_a_in', 'w') as filehandle:
        for entry in contents:
            filehandle.write('%s %s %s\n' % (entry[0], space, entry[1]))


def write_single_site_pot_input(symbol, crystal_type, a, nspins, moment, xc,
                                lmax, print_level, ncomp, conc, mt_radius,
                                ws_radius, egrid, ef, niter, mp):
    """
    Function to write input file for generating single site
    potential using 'newss' command
    Parameters
    ----------

    symbol: str
        Chemical symbol of the element
    crystal_type: int
                1 for FCC, 2 for BCC.
    a: float
        The lattice constant.
    nspins: int
        number of spins.
    moment: float
        Magnetic moment. If nspins = 2 and moment = 0 during input,
        moment will be changed to values from
        this dictionary: {Fe': 2.1, 'Co': 1.4, 'Ni': 0.6}
    xc: int
        ex-cor type (1=vb-hedin,2=vosko).
    lmax: int
        angular momentum quantum number cutoff value.
    print_level: int
        Print level.
    ncomp: int
        Number of components.
    conc: float
        Concentrations.
    mt_radius: float
            mt_radius.
    ws_radius: float
        ws_radius.
    egrid: vector
        e-grid vector of form (ndiv(=#div/0.1Ryd), bott, eimag).
    ef: float
        Estomate of fermi energy.
    niter: int
        Maximum number of SCF iterations.
    mp: float
        Mixing parameter for SCF iterations.
    """

    title = symbol + ' Single Site Potential'
    output_file = symbol + '_ss_out'
    input_file = symbol + '_a_pot'
    pot_file = symbol + '_ss_pot'
    keep_file = symbol + '_ss_k'

    z = atomic_numbers[symbol]
    a = a / Bohr

    if moment == 0.:
        if nspins == 2:
            if symbol in ['Fe', 'Co', 'Ni']:
                moment = magmoms[symbol]

    space = '                    '
    contents = [('Title', title),
                (output_file,
                 'Output file name. If blank, data will show on screen'),
                (print_level, 'Print level'),
                (crystal_type, 'Crystal type (1=FCC,2=BCC)'),
                (lmax, 'lmax'),
                (round(a, 3), 'Lattice constant'),
                (nspins, 'Number of spins'),
                (xc, 'Exchange Correlation type (1=vb-hedin,2=vosko)'),
                (ncomp, 'Number of components'),
                (str(z) + '  ' + str(moment),
                 'Atomic number, Magnetic moment'),
                (conc, 'Concentrations'),
                (str(mt_radius / Bohr) + '  ' + str(ws_radius / Bohr),
                 'mt radius, ws radius'),
                (input_file, 'Input potential file'),
                (pot_file, 'Output potential file'),
                (keep_file, 'Keep file'),
                (str(egrid[0]) + ' ' + str(egrid[1]) + ' ' + str(egrid[2]),
                 'e-grid: ndiv(=#div/0.1Ryd), bott, eimag'),
                (str(round(ef / Rydberg, 3)) + ' '
                 + str(round(ef / Rydberg, 3)),
                 'Fermi energy (estimate)'),
                (str(niter) + ' ' + str(mp),
                 'Number of scf iterations, Mixing parameter')]

    with open(str(symbol) + '_ss_in', 'w') as filehandle:
        for entry in contents:
            filehandle.write('%s %s %s\n' % (entry[0], space, entry[1]))


def write_input_parameters_file(atoms, parameters):
    """Write the main input file for 'mst2' command. This file contains all
    essential input parameters required for calculation"""
    energy_params = ['etol', 'ptol', 'ftol',
                     'offset_energy_pt',
                     'em_switch']  # Parameters with units of energy
    spatial_params = ['liz_cutoff', 'max_core_radius',
                      'max_mt_radius', 'core_radius',
                      'mt_radius']  # Parameters with units of length
    vector_params = ['uniform_grid', 'grid_origin', 'grid_1',
                     'grid_2', 'grid_3', 'grid_pts', 'kpts',
                     'moment_direction', 'constrain_field',
                     'liz_shell_lmax', 'em_mix_param']  # vector parameters
    # Header
    hline = 80 * '='
    separator = 18 * ' ' + 3 * ('* * *' + 14 * ' ')
    header = [hline, '{:^80s}'.format('Input Parameter Data File'),
              hline, separator, hline,
              '{:^80}'.format('System Related Parameters'), hline]

    natoms = ['No. Atoms in System (> 0)  ::  '
              + str(len(atoms)), hline, separator, hline]

    # Get number of atoms from CPA sites if self.parameters['method'] == 3
    if 'method' in parameters.keys():
        if parameters['method'] == 3:
            natoms = ['No. Atoms in System (> 0)  ::  '
                      + str(len(atoms.info['CPA'])), hline, separator, hline]

    header += natoms

    with open('i_new', 'w') as filehandle:
        for entry in header:
            filehandle.write('%s\n' % entry)

    # Rest of the parameters:
    contents = []

    for key in parameters.keys():
        if key in energy_params:
            parameters[key] = parameters[key] / Rydberg

        if key in spatial_params:
            parameters[key] = parameters[key] / Bohr

        if key in vector_params:
            parameters[key] = str(parameters[key])[1:-1]

        if key == 'in_pot':
            for index in parameters['in_pot'].keys():
                contents.append(defaults[key] + '  ::  '
                                + index + ' ' + parameters['in_pot'][index])
        else:
            contents.append(defaults[key] + '  ::  ' + str(parameters[key]))

    with open('i_new', 'a') as filehandle:
        for entry in contents:
            filehandle.write('%s\n' % entry)


def generate_starting_potentials(atoms, crystal_type, a, cpa=False, nspins=1,
                                 moment=0., xc=1, lmax=3, print_level=1,
                                 ncomp=1, conc=1., mt_radius=0., ws_radius=0,
                                 egrid=(10, -0.4, 0.3), ef=9.5, niter=50,
                                 mp=0.1):
    """
    Function to generate single site starting potentials for
    all elements in an atoms object

    Parameters
    ----------
    atoms: The Atoms object to generate starting potentials.
    crystal_type: int
        1 for FCC, 2 for BCC.
    a: float
        The lattice constant.
    cpa: bool
        If True, use atoms.info['CPA'] to generate starting potentials
        for all elements in CPA formalism
    nspins: int
        number of spins.
    moment: float
        Magnetic moment. If nspins = 2 and moment = 0 during input,
        moment will be changed to values from this
        dictionary: {Fe': 2.1, 'Co': 1.4, 'Ni': 0.6}
    xc: int
        ex-cor type (1=vb-hedin,2=vosko).
    lmax: int
        angular momentum quantum number cutoff value.
    print_level: int
        Print level.
    ncomp: int
        Number of components.
    conc: float
        Concentrations.
    mt_radius: float
        mt_radius.
    ws_radius: float
        ws_radius.
    egrid: vector
        e-grid vector of form (ndiv(=#div/0.1Ryd), bott, eimag).
    ef: float
        Estimate of fermi energy.
    niter: int
        Maximum number of SCF iterations.
    mp: float
        Mixing parameter for SCF iterations.
    """
    if cpa:
        species = []
        for site in atoms.info['CPA']:
            for element in site:
                if element == 'index':
                    pass
                else:
                    species.append(element)
        species = np.unique(species)
    else:
        species = np.unique(atoms.get_chemical_symbols())

    for symbol in species:
        # Generate atomic potential
        write_atomic_pot_input(symbol, nspins=nspins, moment=moment,
                               xc=xc, niter=niter, mp=mp)

        subprocess.run('newa', stdin=open(symbol + '_a_in'), check=True)

        # Generate single site potential
        write_single_site_pot_input(symbol=symbol, crystal_type=crystal_type,
                                    a=a, nspins=nspins, moment=moment, xc=xc,
                                    lmax=lmax, print_level=print_level,
                                    ncomp=ncomp, conc=conc,
                                    mt_radius=mt_radius,
                                    ws_radius=ws_radius,
                                    egrid=egrid, ef=ef, niter=niter, mp=mp)

        subprocess.run('newss', stdin=open(symbol + '_ss_in'), check=True)


class MuST(FileIOCalculator):
    """
    Multiple Scattering Theory based ab-initio calculator.
    Capable of performing LSMS, KKR and KKR-CPA calculations.
    """

    implemented_properties = ['energy']
    command = 'mst2 < i_new'

    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label='mst', atoms=None, **kwargs):
        FileIOCalculator.__init__(self, restart, ignore_bad_restart_file,
                                  label, atoms, **kwargs)

        if 'ntasks' in self.parameters:
            self.command = 'mpirun -np ' + str(self.parameters['ntasks']) \
                           + ' ' + self.command

    def clean(self):
        """
        Delete the input parameters file, positions data file
        and default output files generated by MuST
        """

        files = ['o_n00000_D', 'k_n00000_D', 'position.dat', 'i_new']
        for f in files:
            try:
                os.remove(f)
            except OSError:
                pass

    def write_input(self, atoms, properties=None, system_changes=None):
        """ Write input files"""

        FileIOCalculator.write_input(self, atoms,
                                     properties=None, system_changes=None)
        # Write positions using CPA sites if self.parameters['method'] == 3
        if 'method' in self.parameters.keys():
            if self.parameters['method'] == 3:
                method = self.parameters['method']
            else:
                method = None

        else:
            method = None

        write_positions_input(atoms, method=method)
        write_input_parameters_file(atoms=atoms, parameters=self.parameters)

    def read_results(self):
        """ Read results from output files"""
        outfile = glob.glob('k_n00000_D')[0]
        with open(outfile, 'r') as file:
            lines = file.readlines()

        e_offset = float(lines[7].split()[-1])

        results = {tag: value for tag, value in zip(lines[9].split(),
                                                    lines[-1].split())}
        read_energy = (float(results['Energy']) + e_offset)

        convergence = False

        outfile = glob.glob('o_n00000_D')[0]
        with open(outfile, 'r') as file:
            lines = file.readlines()

        for line in lines:
            if 'SCF Convergence is reached' in line:
                convergence = True
                break

        if not convergence:
            raise SCFError()
        else:
            self.results['energy'] = read_energy * Rydberg
