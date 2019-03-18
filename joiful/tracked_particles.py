## File with functions pertaining to creating particle trajectories
import numpy as np
import scipy as sp
import h5py

## The name used for the dataset which stores the timesteps
timeStep_str="timeSteps"

coordinateDict={4:['x','px','py','pz'],
                5:['x','y','px','py','pz'],
                6:['x','y','z','px','py','pz']}


def get_IDs_in_box(trackedDiagnostic, timeStep, box, **kwargs):
    # Function for finding the particles which are in a given
    # phase-space box at a given time step.
    #

    ## Checks the keys
    for key in kwargs.keys():
        if key != '': ## Dummy check, to see that there are no keywords given
            raise KeyError('Key: %s not supported'%key)
    #end key check
    
    particles = trackedDiagnostic.getData(timestep=timeStep)
    ## The size of the array of partilces (in the first dimension of the box) 
    N_particles=particles[timeStep][list(box.keys())[0]].size
    ## Array of bools with the particles satisfying the box conditions
    inBox=np.full(N_particles, True, dtype=bool)

    ## looping over the box dimensions
    for dim in box.keys():
        tmp_coord=particles[timeStep][dim]
        condition=np.sort(box[dim])
        ##Check that only two limits are given for the box dimension
        if condition.size != 2:
            print("Recived box condition %s with %d bounds"%(str(dim),condition.size))
            raise ValueError("The box condition can only have 2 and only 2 bounds!")
        ## Refining the boolean array with the current box condition
        inBox=np.logical_and(inBox,## Particles satisfying the previous criterias, and
                                   ## the ones satisfying the present criteria
                             np.logical_and(tmp_coord>condition[0],tmp_coord<condition[1]))
    #end for keys in box

    if not np.any(inBox):
        print("Did not find any particles in the box. Maybe the conditions are too tight.")
    
    return particles[timeStep]['Id'][inBox]


def extract_trajectories(trackedDiagnostic, IDs, **kwargs):
    # Function for exctracting the particle trajectories of the
    # particles with ID in IDs
    # 
    # trackedDiagnostic - SmileiSimulation.TrackParticles(species=<species>, sort=False)
    #
    # IDs - list([<particle ID's of the particles to be tracked>])

    ## Checks the keys
    for key in kwargs.keys():
        if key=='path':
            path=kwargs[key]
        # elif key=='nbrOutputDims':
        #     nbrOutputDims=kwargs[key]
        elif key=='coordinates':
            coordinates=kwargs[key]
        elif key=='append':
            append=kwargs[key]
        else:
            raise KeyError('Key: %s not supported'%key)
    ## Default values
    if 'path' not in locals():
        path='/tmp'
    ## Then sets the filename
    if path.endswith('.hdf5'):
        filename=path
    else:
        filename=path+'/trajectories_'+trackedDiagnostic.species+'.hdf5'

    if 'coordinates' not in locals():
        coordinates=coordinateDict[4]#['x','px','py','pz']
    nbrOutputDims=len(coordinates)

    ## The time steps
    timeSteps=trackedDiagnostic.getTimesteps()
    N_timeSteps=timeSteps.size

    ## Opens a hdf5 file to write the trajectories to
    if append:
        hf=h5py.File(filename,'a')
    else:
        hf=h5py.File(filename,'w')
    
    ## Creates a dataset with all the timesteps
    if timeStep_str not in list(hf.keys()):
        hf.create_dataset(timeStep_str,data=timeSteps)
    
    ## Creates datasets with the particles
    datasets={}
    new_IDs=[]
    for Id in IDs:
        ## Checks that the particle Id is not already in the file.
        if str(Id) not in list(hf.keys()):
            new_IDs+=[Id]
            datasets[Id]=hf.create_dataset(str(Id),(N_timeSteps,nbrOutputDims),dtype='float64')
    ##end for IDs
    IDs=new_IDs
    ## Checks if the list of IDs is empty, if so return.
    if len(IDs) == 0:
        print("I have nothing to add.")
        return
    
    ## Defining how often to print progress
    printout_every=N_timeSteps//100
    
    for i in range(N_timeSteps):
        n=timeSteps[i]
        ## Printout every 1% of the progress
        if i % printout_every == 0:
            print("i = %4d of %d (%0.0f %% done)" % (i,N_timeSteps,100*i/N_timeSteps))
        ## Opening the file
        particles=trackedDiagnostic.getData(timestep=n)
        for Id in IDs:
            ## Finds the correct particle
            index=np.where(particles[n]['Id']==Id)[0]
            ## Checks that the index is valid
            if index.size>0:
                index=index[0]
            else:
                break
            j=0
            for coord in coordinates:
                datasets[Id][i,j]=particles[n][coord][index]
                j+=1
            ##end for coordinates
        ##end for IDs
    ##end for timesteps

    print("Writing to file")
    hf.close()
    print("DONE!!!")
    return





def read_trajectories(filename, **kwargs):
    for key in kwargs.keys():
        if key=='dataset_IDs':
            dataset_IDs=kwargs[key]
        else:
            raise KeyError('Key: %s not supported'%key)
    ##end for kwargs.keys
    h5File=h5py.File(filename,'r')
    h5Keys=h5File.keys()
    ## Default values
    if 'dataset_IDs' not in locals():
        dataset_IDs=list(h5Keys)
    else:
        ## Converts the dataset_IDs to list of strings
        dataset_IDs=list(dataset_IDs)
        for i in range(len(dataset_IDs)):
            dataset_IDs[i]=str(dataset_IDs[i])
        ## Adds the timesteps to the retrived dataset_IDs if missing
        if timeStep_str not in dataset_IDs:
            dataset_IDs=[timeStep_str]+list(dataset_IDs)
        ## Check that the supplied dataset_IDs are valid
        for c in dataset_IDs:
            if c not in h5Keys:
                print("Unsupported coordinate '%s'."%c)
                print("This file only supports the dataset_IDs:")
                print(list(h5Keys))
                return -1
    ##end default values

    ## Retrive the desired data from file
    dataDict={}
    for Id in dataset_IDs:
        if Id == timeStep_str:
            dataDict[Id] = np.array(h5File.get(Id))
        else:
            data=np.array(h5File.get(Id))
            Ncoord=data.shape[1]
            dataDict[int(Id)]={}
            for i in range(Ncoord):
                dataDict[int(Id)][coordinateDict[Ncoord][i]]=data[:,i]
    ## end retrive data from file
    return dataDict




