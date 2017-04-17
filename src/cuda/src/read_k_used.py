import numpy as np 

def len_kkfile(datafile):
        #Calculates the length of the time series file
        kkac = open(datafile, "r")

        line = kkac.readline() # Skip the title line

        nstep = 0
        for line in kkac:
                nstep = nstep + 1

        kkac.close()

        nstep = nstep - 1 #remove the extra step

        print nstep

        return nstep

def read_kkfile(datafile):
        nstep = len_kkfile(datafile)

        kkac = open(datafile, "r")

        line = kkac.readline()
        line = line.split()

        headers = line

        var_list = np.zeros((nstep+1, len(headers)))

        kk = 0
        for line in kkac:
                line = line.split()
                for ii in np.arange((len(headers))):
                        var_list[kk, ii] = float(line[ii])
                kk = kk+1

        kkac.close()

        return headers, var_list

class kkdata:
        'Includes all the data from the used k file.'
        def __init__(self, kkfile):
                self.headers, self.var_list  = read_kkfile(kkfile)
                #The order of variable in the table: kk_vec_x kk_vec_y kk_vec_z phi forcing_kk_part_x forcing_kk_part_y forcing_kk_part_z
                for idx, name in enumerate(self.headers, start=0):
                        if name == "kk_vec_x":
                                self.kk_vec_x = self.var_list[:,idx]
                        elif name == "kk_vec_y":
                                self.kk_vec_y = self.var_list[:,idx]
                        elif name == "kk_vec_z":
                                self.kk_vec_z = self.var_list[:,idx]
                        elif name == "phi":
                                self.phi = self.var_list[:,idx]
                        elif name == "forcing_kk_part_x":
                                self.forcing_kk_part_x = self.var_list[:,idx]
                        elif name == "forcing_kk_part_y":
                                self.forcing_kk_part_y = self.var_list[:,idx]
                        elif name == "forcing_kk_part_z":
                                self.forcing_kk_part_z = self.var_list[:,idx]
                        else:
                                print "Error: " + name + " not recognized is a time series header!"

