#import matplotlib
import numpy as np
import time
from math import *

class data_set():
    def __init__(self):
        self.parent = 0
        self.decay_options = np.array([])
        self.type = np.array([])
        self.event_index = np.array([])
        self.particle_index = np.array([])
        self.pid = np.array([])
        self.p = np.array([])
        self.pT = np.array([])
        self.eta = np.array([])
        self.phi = np.array([])
        self.m = np.array([])
        self.tau = np.array([])
        self.xi = np.array([])
        self.xf = np.array([])
        self.xPV = np.array([])
    
    def print_vbls(self):
        print("0:type")
        print("1:event_index")
        print("2:particle_index")
        print("3:pid")
        print("4:p")
        print("5:pT")
        print("6:eta")
        print("7:phi")
        print("8:m")
        print("9:tau")
        print("10:xi")
        print("11:xf")
        print("12:xPV")

    def save(self, prefix):
        #Reset 4-vector sets
        p = np.reshape(self.p, (self.p.shape[0]*4))
        xi = np.reshape(self.xi, (self.xi.shape[0]*4))
        xf = np.reshape(self.xf, (self.xf.shape[0]*4))
        xPV = np.reshape(self.xPV, (self.xPV.shape[0]*4))
        #Append all values
        full_set = np.array([])
        full_set = np.append(full_set, np.array([self.parent, self.decay_options.shape[0]]))
        full_set = np.append(full_set, self.decay_options)
        full_set = np.append(full_set, self.type)
        full_set = np.append(full_set, self.event_index)
        full_set = np.append(full_set, self.particle_index)
        full_set = np.append(full_set, self.pid)
        full_set = np.append(full_set, p)
        full_set = np.append(full_set, self.pT)
        full_set = np.append(full_set, self.eta)
        full_set = np.append(full_set, self.phi)
        full_set = np.append(full_set, self.m)
        full_set = np.append(full_set, self.tau)
        full_set = np.append(full_set, xi)
        full_set = np.append(full_set, xf)
        full_set = np.append(full_set, self.xPV)
        np.save(prefix + "_set.npy", full_set)
        return full_set

    def clear_data(self):
        self.type = np.array([])
        self.event_index = np.array([])
        self.particle_index = np.array([])
        self.pid = np.array([])
        self.p = np.array([])
        self.pT = np.array([])
        self.eta = np.array([])
        self.phi = np.array([])
        self.m = np.array([])
        self.tau = np.array([])
        self.xi = np.array([])
        self.xf = np.array([])
        self.xPV = np.array([])

    def init_zero_set(self):
        self.type = np.array("-")
        self.event_index = np.array([-1.])
        self.particle_index = np.array([-1.])
        self.pid = np.array([0.])
        self.p = np.array([[0.,0.,0.,0.]])
        self.pT = np.array([0.])
        self.eta = np.array([0.])
        self.phi = np.array([0.])
        self.m = np.array([0.])
        self.tau = np.array([0.])
        self.xi = np.array([[0.,0.,0.,0.]])
        self.xf = np.array([[0.,0.,0.,0.]])
        self.xPV = np.array([[0.,0.,0.,0.]])

    def size(self):
        return self.pid.shape[0]

    def set_parent(self, parent_pid):
        self.parent = parent_pid

    def set_decay_options(self, dopts):
        self.decay_options = dopts
    
    def fill_custom(self, fileName):
        start_time = time.time()
        f = open(fileName)

        for line in f.readlines():
            l = line.split();
            self.type = np.append(self.type, int(l[0]))
            self.event_index = np.append(self.event_index, int(l[1]))
            self.particle_index = np.append(self.particle_index, int(l[2]))
            self.pid = np.append(self.pid, int(l[3]))
            self.p = np.append(self.p,np.array([float(l[10]),float(l[4]),float(l[5]),float(l[6])]))
            self.pT = np.append(self.pT, float(l[7]))
            self.eta = np.append(self.eta, float(l[8]))
            phi = float(l[9])
            if phi < 0.: phi = phi + 2*np.pi ### ????
            self.phi = np.append(self.phi, phi)
            self.m = np.append(self.m, float(l[11]))
            self.tau = np.append(self.tau, float(l[12]))
            self.xi = np.append(self.xi, np.array([float(l[13]),float(l[14]),float(l[15]), float(l[16])]))
            self.xf = np.append(self.xf, np.array([float(l[17]),float(l[18]),float(l[19]),float(l[20])]))
            self.xPV = np.append(self.xPV, np.array([float(l[21]),float(l[22]),float(l[23]),float(l[24])]))

        f.close()
        
        self.p = np.reshape(self.p,(int(self.p.shape[0]/4),4))
        self.xi = np.reshape(self.xi,(int(self.xi.shape[0]/4),4))
        self.xf = np.reshape(self.xf,(int(self.xf.shape[0]/4),4))
        self.xPV = np.reshape(self.xPV,(int(self.xPV.shape[0]/4),4))
        print("Fill custom: --- %s seconds ---" % (time.time()-start_time))
    
    def add_HEPMC_particle(self, eid, l, count):
        self.type = np.append(self.type, "N")
        self.event_index = np.append(self.event_index, eid)
        self.particle_index = np.append(self.particle_index, int(l[1]))
        self.pid = np.append(self.pid,int(l[3]))
        px, py, pz, e = float(l[4]), float(l[5]), float(l[6]), float(l[7])
        self.pT = np.append(self.pT, np.sqrt(px*px + py*py))
        
        if(np.sqrt(px*px + py*py) > 0):
            
            self.eta = np.append(self.eta, np.arctanh(pz/(np.sqrt(px*px + py*py + pz*pz)))) 
            phi = np.arcsin(np.abs(py)/np.sqrt(px*px + py*py)) 
            if(px > 0. and py < 0.):#Adding phi in 0-2pi angle measure
                phi = 2*np.pi - phi
            elif(px < 0. and py > 0.):
                phi = np.pi - phi
            elif(px < 0. and py < 0.):
                phi = phi + np.pi
            self.phi = np.append(self.phi, phi)

        else: #Not a problem for Ap decays, just Higgs->Ap
            self.eta = np.append(self.eta, 1000.) #Come back to this 
            self.phi = np.append(self.phi, 0.) #Really just N/A 
        
        self.m = np.append(self.m, np.sqrt(e*e - px*px - py*py - pz*pz)) 
        self.tau = np.append(self.tau, 0.) #Come back to this 
        
        if(count > 1):
            self.p = np.append(self.p, np.array([[e,px,py,pz]]),axis=0) #Have to add element first
            self.xi = np.append(self.xi, np.array([[0.,0.,0.,0.]]),axis=0) #Gets edited when vertex found
            self.xf = np.append(self.xf, np.array([[0.,0.,0.,0.]]),axis=0) #Gets edited when vertex found
            self.xPV = np.append(self.xPV, np.array([[0.,0.,0.,0.]]),axis=0) #Gets edited when vertex found
        
        elif(count == 1):
            self.p[0] = np.array([[e,px,py,pz]])

    def fill_HEPMC(self,fileName):
        start_time = time.time()
        #First, correctly initialize x and p arrays
        #There is probably a neater way of doing this
        self.p = np.array([[0.,0.,0.,0.]])
        self.xi = np.array([[0.,0.,0.,0.]])
        self.xf = np.array([[0.,0.,0.,0.]])
        self.xPV = np.array([[0.,0.,0.,0.]])

        #Open hepmc data file
        f = open(fileName)
        
        #Temporary storage
        temp_eid = -1 #So that it starts at 0.
        vid = np.array([])
       
        FOUND = 0 #Counts up number of Ap per event
        count_particles = 0
        
        for line in f.readlines():
            l = line.split();
         
            if line.startswith("E"):
                temp_eid = temp_eid + 1
                FOUND = 0 #Reset at the beginning of each new event
    
            elif line.startswith("P"):
               
                pid = int(l[3])
                if pid == self.parent:#If it's an Ap
                    FOUND = FOUND + 1
                    vid= np.append(vid,0)#Create an empty vertex location
                    count_particles = count_particles + 1
                    self.add_HEPMC_particle(temp_eid, l, count_particles)
                
                elif np.where(np.isin(self.decay_options, np.abs(pid)))[0].shape[0] > 0: #If it's on the decay options list
                    if FOUND > 0:#If an Ap has already been found (just to save time)
                        v_prod = int(l[2]) #Vertex mother of particle
                        vertex_match = np.where(np.isin(vid, v_prod))[0]
                       
                        if vertex_match.shape[0] > 0 and temp_eid == self.event_index[vertex_match[-1]]:
                            vid= np.append(vid,0)
                            count_particles = count_particles + 1
                            self.add_HEPMC_particle(temp_eid, l, count_particles)
                            self.xi[-1] = self.xf[vertex_match[-1]]
            
            elif line.startswith("V"):
            
                if(len(l) > 4 and l[4] == '@' and l[3].rfind(",") == -1):#If displaced decay found 
                    decayed = int(l[3].strip('[]'))#Get the index of the decay particle
                
                    if FOUND > 0:#Check if current event is associated with an Ap
                        index_match = np.where(np.isin(self.particle_index, decayed))[0]#Locate Aps in array
                        
                        if(index_match.shape[0] > 0):
                            hyp_index = index_match[-1] #Index hypothesis

                            if self.event_index[hyp_index] == temp_eid and self.pid[hyp_index] == self.parent:
                                vid[hyp_index] = int(l[1]) #Fill vid with the last match
                                self.xf[hyp_index] = np.array([float(l[8]),float(l[5]),float(l[6]),float(l[7])])
            
        self.xi = self.xi*0.1
        self.xf = self.xf*0.1
        self.xPV = self.xPV*0.1
        print("Fill HEPMC: --- %s seconds ---" % (time.time()-start_time))
    
    def fill_saved(self, fileName):
        start_time = time.time()
        
        #Load and initialize
        saved = np.load(fileName)
        self.parent = saved[0]
        self.decay_options = saved[2:int(2+saved[1])]
        saved = saved[int(2+saved[1]):]

        #Fill
        nentries = int(saved.shape[0]/25)
        self.type = saved[0:nentries]; saved = saved[nentries:]
        self.event_index = saved[0:nentries]; saved = saved[nentries:]
        self.particle_index = saved[0:nentries]; saved = saved[nentries:]
        self.pid = saved[0:nentries]; saved = saved[nentries:]
        self.p = saved[0:4*nentries]; saved = saved[4*nentries:]
        self.pT = saved[0:nentries]; saved = saved[nentries:]
        self.eta = saved[0:nentries]; saved = saved[nentries:]
        self.phi = saved[0:nentries]; saved = saved[nentries:]
        self.m = saved[0:nentries]; saved = saved[nentries:]
        self.tau = saved[0:nentries]; saved = saved[nentries:]
        self.xi = saved[0:4*nentries]; saved = saved[4*nentries:]
        self.xf = saved[0:4*nentries]; saved = saved[4*nentries:]
        self.xPV = saved[0:4*nentries]; saved = saved[4*nentries:]
        
        #Reshape necessary elements
        self.p = np.reshape(self.p,(int(self.p.shape[0]/4),4))
        self.xi = np.reshape(self.xi,(int(self.xi.shape[0]/4),4))
        self.xf = np.reshape(self.xf,(int(self.xf.shape[0]/4),4))
        self.xPV = np.reshape(self.xPV,(int(self.xPV.shape[0]/4),4))
        print("Fill saved: --- %s seconds ---" % (time.time()-start_time))
        

    def get(self, mask): 
        data = data_set()
        data.type = self.type[mask]
        data.event_index = self.event_index[mask]
        data.particle_index = self.particle_index[mask]
        data.pid = self.pid[mask]
        data.p = self.p[mask]
        data.pT = self.pT[mask]
        data.eta = self.eta[mask]
        data.phi = self.phi[mask]
        data.m = self.m[mask]
        data.tau = self.tau[mask]
        data.xi = self.xi[mask]
        data.xf = self.xf[mask]
        data.xPV = self.xPV[mask]
        return data
    
    def remove(self, mask):
        return self.get(~mask)
   
    def get_type(self, my_type):
        mask = self.type == my_type
        return self.get(mask)
    
    def remove_type(self, my_type):
        mask = self.pid != my_type
        return self.get(mask)

    def get_pid(self, p_id):
        mask = self.pid == p_id
        return self.get(mask)
    
    def remove_pid(self, p_id):
        mask = self.pid != p_id
        return self.get(mask)
    
    def get_anti(self):
        mask = self.pid < 0.
        return self.get(mask)
    
    def remove_anti(self):
        mask = self.pid > 0.
        return self.get(mask)

    def copy(self):
        mask = np.ones(self.size(),dtype=bool)
        return self.get(mask)

    def reset_elements(self, mask):
        #Try to think about ways of combining this with set
        self.event_index[mask] = 0.
        self.particle_index[mask] = 0.
        self.pid[mask] = 0.
        self.p[mask] = 0.
        self.pT[mask] = 0.
        self.eta[mask] = 0. 
        self.phi[mask] = 0.
        self.m[mask] = 0.
        self.tau[mask] = 0.
        self.xi[mask] = 0.
        self.xf[mask] = 0.
        self.xPV[mask] = 0.
        return self
    
    def set_elements(self, mask, data):
        self.event_index[mask] = data.event_index
        self.particle_index[mask] = data.particle_index
        self.pid[mask] = data.pid
        self.p[mask] = data.p
        self.pT[mask] = data.pT
        self.eta[mask] = data.eta 
        self.phi[mask] = data.phi
        self.m[mask] = data.m
        self.tau[mask] = data.tau
        self.xi[mask] = data.xi
        self.xf[mask] = data.xf
        self.xPV[mask] = data.xPV
        return self
    
    def delete(self, i_elem):
        #Try to think about ways of combining this with remove
        self.type = np.delete(self.type, i_elem)
        self.event_index = np.delete(self.event_index, i_elem)
        self.particle_index = np.delete(self.particle_index, i_elem)
        self.pid = np.delete(self.pid, i_elem)
        self.p = np.delete(self.p, i_elem, axis=0)
        self.pT = np.delete(self.pT, i_elem)
        self.eta = np.delete(self.eta, i_elem)
        self.phi = np.delete(self.phi, i_elem)
        self.m = np.delete(self.m, i_elem)
        self.tau = np.delete(self.tau, i_elem)
        self.xi = np.delete(self.xi, i_elem, axis=0)
        self.xf = np.delete(self.xf, i_elem, axis=0)
        self.xPV = np.delete(self.xPV, i_elem, axis=0)
        return self
    
    def append(self, new_data):
        self.type = np.append(self.type, new_data.type)
        self.event_index = np.append(self.event_index, new_data.event_index)
        self.particle_index = np.append(self.particle_index, new_data.particle_index)
        self.pid = np.append(self.pid, new_data.pid)
        self.p = np.append(self.p, new_data.p, axis=0)
        self.pT = np.append(self.pT, new_data.pT)
        self.eta = np.append(self.eta, new_data.eta)
        self.phi = np.append(self.phi, new_data.phi)
        self.m = np.append(self.m, new_data.m)
        self.tau = np.append(self.tau, new_data.tau)
        self.xi = np.append(self.xi, new_data.xi, axis=0)
        self.xf = np.append(self.xf, new_data.xf, axis=0)
        self.xPV = np.append(self.xPV, new_data.xPV, axis=0)
        return self
    
    def copy_single(self, vbl, size, i):
        #Unfinished
        new_set = data_set()
        
        if(vbl == 0):
            new_set.type = np.empty(size); new_set.type.fill(self.type[i])
        elif(vbl == 1):
            new_set.type = np.empty(size); new_set.type.fill(self.type[i])
        elif(vbl == 2):
            new_set.type = np.empty(size); new_set.type.fill(self.type[i])
        elif(vbl == 3):
            new_set.type = np.empty(size); new_set.type.fill(self.type[i])
        elif(vbl == 4):
            #*
            new_set.p = np.reshape(np.zeros(size*4), (size, 4))
            new_set.p[:,0].fill(self.p[i][0])
            new_set.p[:,1].fill(self.p[i][1])
            new_set.p[:,2].fill(self.p[i][2])
            new_set.p[:,3].fill(self.p[i][3])
        elif(vbl == 5):
            new_set.type = np.empty(size); new_set.type.fill(self.type[i])
        elif(vbl == 6):
            #Note: Only eta/phi have been fully entered
            new_set.eta = np.empty(size); new_set.eta.fill(self.eta[i])
        elif(vbl == 7):
            #*
            new_set.phi = np.empty(size); new_set.phi.fill(self.phi[i])
        elif(vbl == 8):
            new_set.type = np.empty(size); new_set.type.fill(self.type[i])
        elif(vbl == 9):
            new_set.type = np.empty(size); new_set.type.fill(self.type[i])
        elif(vbl == 10):
            new_set.type = np.empty(size); new_set.type.fill(self.type[i])
        elif(vbl == 11):
            new_set.type = np.empty(size); new_set.type.fill(self.type[i])
        elif(vbl == 12):
            new_set.type = np.empty(size); new_set.type.fill(self.type[i])

        return new_set
    
    def find_repeats(self):
        #Find all pairs of repeated events
        repeats = np.zeros(self.size(),dtype=bool)

        #for i in range(self.size()-1):
         #   if self.event_index[i] == self.event_index[i+1]:
          #      repeats[i] = True
           #     repeats[i+1] = True
            #    i = i+2
        i = 0
        f = self.size()
        while (i < f):
            if self.event_index[i] == self.event_index[i+1]:
                repeats[i] = True
                repeats[i+1] = True
                i = i+2  
            else: 
                i = i+1
                if(i == self.size()-1): i = i+1 #Method for exiting loop

        return repeats
    
    def create_even_odd_split(self):
        #Split along even and odd entries
        split_mask = np.ones(self.size(),dtype=bool)
        
        for i in range(self.size()):
            if i%2 == 1:
                split_mask[i] = False

        return split_mask

    def split_set(self, split_mask):
        split1 = self.get(split_mask)
        split2 = self.get(~split_mask)
        return split1, split2

    def is_equal_to(self, data):
        mask = (self.type != data.type)
        if(mask[mask].shape[0] > 0): return False
        
        mask = self.event_index != data.event_index
        if(mask[mask].shape[0] > 0): return False
        
        mask = self.particle_index != data.particle_index
        if(mask[mask].shape[0] > 0): return False
        
        mask = self.pid != data.pid
        if(mask[mask].shape[0] > 0): return False
        
        mask = self.p != data.p
        if(mask[mask].shape[0] > 0): return False
        
        mask = self.pT != data.pT
        if(mask[mask].shape[0] > 0): return False
        
        mask = self.eta != data.eta
        if(mask[mask].shape[0] > 0): return False
        
        mask = self.phi != data.phi
        if(mask[mask].shape[0] > 0): return False
        
        mask = self.m != data.m
        if(mask[mask].shape[0] > 0): return False
        
        mask = self.tau != data.tau
        if(mask[mask].shape[0] > 0): return False
        
        mask = self.xi != data.xi
        if(mask[mask].shape[0] > 0): return False
        
        mask = self.xf != data.xf
        if(mask[mask].shape[0] > 0): return False
        
        mask = self.xPV != data.xPV
        if(mask[mask].shape[0] > 0): return False

        return True












