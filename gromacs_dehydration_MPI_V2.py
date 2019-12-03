import sys
import os
import shutil
import subprocess
import random

class Dehydration:

    #newTempSeries = TempSeries(start_pdb, itp, top, newfilename, moleculetypes, numberofthreads, temperature, waters_remove, segment_time, number_of_seqments)

    def __init__(self, start_pdb, itp_file, top_file, new_filename, moleculetypes, numberofthreads, temperature, waters_remove, segment_time, number_of_seqments, barostat, savedata):
        self._new_filename = new_filename
        self._start_pdb =  start_pdb
        self._itp_file = itp_file
        self._top_file = top_file
        self._moleculetypes = moleculetypes
        self._moleculetypesstring = ''
        for mol in self._moleculetypes:
            self._moleculetypesstring += ' '+mol
        self._numberofthreads = numberofthreads
        self._temperature = temperature
        self._waters_remove = waters_remove
        self._segment_time = segment_time
        self._number_of_seqments = number_of_seqments
        self._barostat = barostat
        self._savedata = savedata

    def run_process(self,exe):
        p = subprocess.Popen(exe, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,bufsize = 10000000, shell=True)
        while(True):
            retcode = p.poll() #returns None while subprocess is running
            line = p.stdout.readline()
            yield line
            if(retcode is not None):
                break

    def run_gmx(self,exe):
        output = ''
        for line in self.run_process(exe):
            #output.append(line)
            output += line             
        return output

    def NPT_write_mdp(self, filename, temperature, picoseconds, generate_velocities, restart, barostat):

        out = ""
        out += "; GAFF NPT dehydration \n"
        out += "; Run parameters\n"
        out += "integrator      = md            ; leap-frog integrator\n"

        #out += "nsteps          = 50000         ; 1 * 50000 = 50 ps\n"
        out += "nsteps          =" + str(picoseconds*1000) + "; steps\n"

        out += "dt              = 0.001         ; 1 fs\n"
        out += "; Output control\n"
        out += "nstxout         = " + str(self._savedata) + "          ; save coordinates every 1.0 ps\n"
        out += "nstvout         = " + str(self._savedata) + "          ; save velocities every 1.0 ps\n"
        out += "nstenergy       = " + str(self._savedata) + "          ; save energies every 1.0 ps\n"
        out += "nstlog          = " + str(self._savedata) + "          ; update log file every 1.0 ps\n"
        out += "; Bond parameters\n"

        #out += "continuation    = no            ; first dynamics run\n"
        out += "continuation    =" + restart + "      ; first dynamics run\n"

        out += "lincs_iter      = 1             ; accuracy of LINCS\n"
        out += "lincs_order     = 4             ; also related to accuracy\n"
        out += "; Neighborsearching\n"
        out += "cutoff-scheme   = Verlet\n"
        out += "ns_type         = grid          ; search neighboring grid cells\n"
        out += "nstlist         = 10            ; 20 fs, largely irrelevant with Verlet\n"
        out += "rcoulomb        = 1.0           ; short-range electrostatic cutoff (in nm)\n"
        out += "rvdw            = 1.0           ; short-range van der Waals cutoff (in nm)\n"
        out += "; Electrostatics\n"
        out += "coulombtype     = PME           ; Particle Mesh Ewald for long-range electrostatics\n"
        out += "vdwtype         = PME           ; Particle Mesh Ewald for long-range van der waals\n"
        out += "pme_order       = 4             ; cubic interpolation\n"
        out += "fourierspacing  = 0.16          ; grid spacing for FFT\n"
        out += "; Temperature coupling is on\n"
        if barostat == 'berendsen':
            out += "tcoupl          = V-rescale   ;  thermostat\n"
        elif barostat == 'Parrinello-Rahman':
            out += "tcoupl          = nose-hoover  ;  thermostat\n"
        out += "tc-grps         =" + self._moleculetypesstring + " ; two coupling groups - more accurate\n"
        #out += "tau_t           = 0.1 0.1       ; time constant, in ps\n"
        tau_t = "tau_t           ="
        for item in self._moleculetypes:
            tau_t += " 0.2"
        tau_t += "; time constant, in ps\n"
        out += tau_t

        #out += "ref_t           = 50  50        ; reference temperature, one for each group, in K\n"
        ref_t = "ref_t           ="
        for item in self._moleculetypes:
            ref_t += " "+str(temperature)
        ref_t += "; reference temperature, one for each group, in K\n"
        out += ref_t

        out += "; Pressure coupling is on\n"
        out += "pcoupl          = " + barostat + " ; Barostat\n"
        out += "pcoupltype      = anisotropic       ; anisotropic scaling of box vectors\n"
        if barostat == 'berendsen':
            out += "tau_p               = 2.0 2.0 2.0 2.0 2.0 2.0                  ; time constant, in ps\n"
        elif barostat == 'Parrinello-Rahman':
            out += "tau_p               = 10.0 10.0 10.0 10.0 10.0 10.0                  ; time constant, in ps\n"
        out += "ref_p               = 1.0 1.0 1.0 0.0 0.0 0.0  ; reference pressure, in bar\n"
        out += "compressibility     = 4.5e-5 4.5e-5 4.5e-5 4.5e-5 4.5e-5 4.5e-5 ; isothermal compressibility of water, bar^-1\n"
        out += "refcoord_scaling    = no\n"
        out += "; Periodic boundary conditions\n"
        out += "pbc             = xyz               ; 3-D PBC\n"
        out += "; Dispersion correction\n"
        out += "DispCorr        = EnerPres          ; account for cut-off vdW scheme\n"
        out += "; Velocity generation\n"
        out += "gen_vel         =" + generate_velocities + " ; assign velocities from Maxwell distribution\n"
        if generate_velocities == 'yes':
            out += "gen_temp        = "+str(temperature)  + "            ; temperature for Maxwell distribution\n"
            out += "gen_seed        = -1            ; generate a random seed\n"

        f = open(filename, 'w')
        f.write(out)
        f.close()

    def check_error(self, string):
        # check if there is error from gromacs and then quit

        if 'Error' in string:
            print '############################### ERROR IN GROMACS ERROR                   ###############################'
            print string
            print '############################### FATAL EXIT                               ###############################'
            exit() 


    def remove_water(self, watertoremove, watertype, grofilename):

        removed_water = 0
   
        ### read .gro file 
        grofile = []
        f = open('./dehydration/'+grofilename, 'r')
        for line in f:
            grofile.append(line) 
        f.close()
    
        name = grofile[0]
        del grofile[0]
        atom_nr = grofile[0]
        del grofile[0]
        unit_cell = grofile[-1]
        del grofile[-1]

        ### read .top file
        topfile = []
        f = open('./dehydration/'+self._top_file, 'r')
        for line in f:
            topfile.append(line)
        f.close()

        # get molecules
        all_molecules = []
        molecule = []
    
        mol = grofile[0].split()[0]
        for line in grofile:
            newline = line.split()
            if newline[0] != mol:
                all_molecules.append(molecule)
                molecule = []
                molecule.append(line)
                mol = line.split()[0]
            else:
                molecule.append(line)
        all_molecules.append(molecule) 

        # sort molecules
        other_molecules = []
        water_molecules = []
        for molecule in all_molecules:
            water = False
            for atom in molecule:
                if watertype in atom:
                    water = True
            if water == True:
                water_molecules.append(molecule)
            else:
                other_molecules.append(molecule)
        print len(other_molecules), len(water_molecules)

        # delete water
#        water_molecules = random.shuffle(water_molecules)

        for i in range(0,watertoremove):
            if len(water_molecules) != 0 and len(water_molecules) > 0:
                del_water = random.randint(0,len(water_molecules)-1)
                print 'delete water:',del_water+1,len(water_molecules)
                del water_molecules[del_water]
                removed_water += 1
        
        # make new .gro file
        all_molecules = []      
        all_molecules_temp = other_molecules + water_molecules
        for moleculetype in self._moleculetypes:
            for molecule in all_molecules_temp:
                if moleculetype in molecule[0]:
                    all_molecules.append(molecule)
            
        atom_number = 1
        mol_number = 1

        new_gro = []

        for molecule in all_molecules:
            for atom in molecule:

                MOL_NUMBER = atom[0:5]
                MOL_TYPE = atom[5:10]
                ATOM_TYPE = atom[10:15]
                ATOM_NUMBER = atom[10:20]

                XYZ_VELOCITY = atom[20:-1] 
                
                MOL_NUMBER = str(mol_number)
                for i in range(0,5-len(MOL_NUMBER)):
                    MOL_NUMBER = ' ' + MOL_NUMBER

                ATOM_NUMBER = str(atom_number)
                for i in range(0,5-len(ATOM_NUMBER)):
                    ATOM_NUMBER = ' ' + ATOM_NUMBER

                NEWATOM = MOL_NUMBER + MOL_TYPE + ATOM_TYPE + ATOM_NUMBER + XYZ_VELOCITY + '\n'
                new_gro.append(NEWATOM)

                atom_number += 1

            mol_number += 1                

        atom_nr = str(int(atom_nr) - 3*removed_water) + '\n'
        new_gro = [name] + [atom_nr] + new_gro + [unit_cell]

        # write new .gro file
        f = open( './dehydration/'+grofilename, 'w')        
        for line in new_gro:
            f.write(line)
        f.close()

        # edit .top file
        newtop = []
        for line in topfile:
            if watertype in line:
                line = line.split()
                number = int(line[1]) - removed_water
                if number > 0:
                    newline = line[0]+'    '+str(number)+'\n'
                    newtop.append(newline)
                else: 
                    # fix molecule types for .mdp file if there is no more water left
                    if watertype in self._moleculetypes: self._moleculetypes.remove(watertype)
                    self._moleculetypesstring = ''
                    for mol in self._moleculetypes:
                        self._moleculetypesstring += ' '+mol        
            else:
                newtop.append(line)

        f = open('./dehydration/'+self._top_file, 'w')
        for line in newtop:
            f.write(line)
        f.close()


    def run(self):
        # run dehydration

        print '############################### Start GROMACS Crystal dehydartion script ###############################'
        print '############################### Made by Anders S. Larsen                 ###############################'

        # setup dirs
        os.mkdir("./dehydration")
        #shutil.copyfile("./energy_minimization/"+energy_minstructure,"./tempseries/"+energy_minstructure)

        shutil.copyfile("./"+self._start_pdb,"./dehydration/"+self._start_pdb)
        shutil.copyfile("./"+self._itp_file,"./dehydration/"+self._itp_file)
        shutil.copyfile("./"+self._top_file,"./dehydration/"+self._top_file)

        # eq loop
        first = True
        laststructure = "NONE" # not used in first loop


        for i in range(0,self._number_of_seqments):

            print '############################### Dehydration segment '+str(i)+'                ###############################'

            # remove water
            if first == True:
                pass
            else:
                self.remove_water(self._waters_remove, 'SOL', laststructure)
                os.remove('./mdout.mdp')    
        
            NPTname = self._new_filename + '_' + str(i) + '.mdp'
            NPToutput = NPTname.replace('.mdp','.gro')
            NPTtpr = NPTname.replace('.mdp','.tpr')

            # make input
            self.NPT_write_mdp("./dehydration/dehydration.mdp", self._temperature, self._segment_time, 'no', 'yes', self._barostat)  

            if first == True:
                call = "gmx_mpi grompp -f ./dehydration/" + 'dehydration.mdp' + " -c ./dehydration/" + self._start_pdb + " -p ./dehydration/" + self._top_file + " -o ./dehydration/"+ NPTtpr
            else:
                call = "gmx_mpi grompp -f ./dehydration/" + 'dehydration.mdp' + " -c ./dehydration/" + laststructure + " -p ./dehydration/" + self._top_file + " -o ./dehydration/"+ NPTtpr
            print 'DO',call
            output = self.run_gmx(call)
            self.check_error(output)
            print output

            # run md 
            #call = "gmx mdrun -nt " + str(self._numberofthreads) + " -deffnm ./dehydration/" + NPTtpr.replace('.tpr','')
            call = "mpirun -np " + str(self._numberofthreads) + " -perhost 16 gmx_mpi mdrun -deffnm ./dehydration/" + NPTtpr.replace('.tpr','')
            print 'DO',call

            output = self.run_gmx(call)
            self.check_error(output)
            print output

            #os.remove('./mdout.mdp')    
            
            laststructure = NPToutput
            first = False

start_pdb = 'THEOPH02_production300K.gro' # can be .gro or .pdb
itp = 'theophylline_GAFF2_RESP.itp'
top = 'theophylline_GAFF2_RESP.top'
moleculetypes = ['MOL','SOL'] # defined in .top file
barostat = 'berendsen'
temperature = 375 # temperature in kelvin
waters_remove = 3 # number of water molecules to remove at once
segment_time = 100 # time between removing water in ps
number_of_seqments = 360
savedata = 5000 # save coordinates every example savedata = 5000 is save data every 5.0 ps
numberofthreads = 32 # number of cpu/threads
newfilename = "THEOPH02" # name for files produced with gromacs

# start simulation
newdehydration = Dehydration(start_pdb, itp, top, newfilename, moleculetypes, numberofthreads, temperature, waters_remove, segment_time, number_of_seqments, barostat, savedata)
newdehydration.run()



