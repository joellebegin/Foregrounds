import numpy as np 
from astropy.cosmology import WMAP9 as cosmo
from astropy.coordinates import Distance
import astropy.units as u
from astropy.cosmology import z_at_value

class Foregrounds:
    ''' 
    generates 3d box with extragalactic foregrounds. Each pixel has a
    specified number of foregrounds, where every foreground has a different 
    spectrum T(nu) = A*v**(-alpha), where A is as specified in generate_amplitudes().

    Parameters
    ----------
    n: int
        dimensions of box. by default, a box 3D box with n**3 pixels is produced
    
    foregrounds_per_pixel: int
        number of foregrounds in one pixel
    
    foreground_alpha: array or tuple of ints or floats
        specifies the gaussian distribution from which the index of the foreground 
        is drawn. Array/tuple with [mean, std]

    mean_150MHz: float 
        Mean temperature at 150MHz  
        

    luminosty_alpha: float 
        index describing the power law distribution from which foreground
        luminosities will be drawn. Very important that it is float and not int, 
        the pareto function does not like ints. 

    
    Methods
    -------
    foreground.frequency_space(bandwidth):
        after having instanciated a foreground object, frequency_space gives a
        uniformly spaced grid in frequency with foreground contamination

    foreground.real_space(central_frequency):
        after having instanciated a foreground object, given some central frequency,
        returns a uniformly spaced grid in real space centered on the distance 
        corresponding to the central frequency.
    '''

    def __init__(self, n=200, foregrounds_per_pix=5, foreground_alpha = [2.5,0.5], 
                mean_150Mhz = 300., luminosity_alpha = 2.75):
        
        self.n = n 
        self.n_f = foregrounds_per_pix
        self.f_alpha = foreground_alpha
        self.mean_150Mhz = mean_150Mhz
        self.l_alpha = luminosity_alpha
        self.NU_21CM = 1420 #MHz


    #========================METHODS THAT OUTPUT BOXES========================#

    def frequency_space(self, bandwith = [100,200]):
        '''
        Returns a box that is uniformly spaced in frequency. The range of frequencies
        is specified by the bandwidth, and the spacing by the size of the box.
        Every slice contains the foreground temperature at a specific frequency. 
        That is, foreground_temps[i] is what the sky would look like at 
        a frequency of self.nu[i]. 

        Parameters
        ----------
        -bandwith: Array/tuple of ints or floats
            range of frequencies that will be considered. 
            Array/tuple with [nu_min, nu_max]
            
        Returns
        -------
        -foreground_temps: 3D numpy array
            3D box where every slice contains foreground temperatures at different
            frequencies
        -nu: 1D numpy array
            The frequencies corresponding to each slice. foreground_temps[i]
            is at frequency nu[i]
        '''

        self.bandwith = bandwith

        #foregrounds in frequency space
        self.generate_frequency_array()
        self.generate_foregrounds(self.nu_arr)

        return (self.foreground_temps, self.nu)

    def real_space(self, central_freq = 150, L = 300):
        '''
        Returns a box that is uniformly spaced in comoving distance. The grid 
        spacing is specified by L and n. 

        Parameters
        ----------
        -central_freq: int or float
            the frequency that will be placed at the center of the box in MHz
        -L: int or float. 
            real space length of the box in Mpc
            
        Returns
        -------
        -foreground_temps: 3D numpy array
            3D box where every slice contains foreground temperatures at different
            frequencies
        '''
        
        self.L = L
        
        #finds distance that central frequency corresponds to
        self.r_0 = self.z_to_dist(self.freq_to_z(central_freq))
        
        self.generate_distances() #finds the comoving distances are contained in box
        self.distances_to_freq() #converts distances back to frequencies
        self.generate_foregrounds(self.frequency_grid)
        
        return self.foreground_temps


    #=============FUNCTIONS RELATED TO THE FREQUENCY SPACE OUTPUT=============#

    def generate_frequency_array(self):
        '''generates the range of frequencies considered.'''

        self.nu = np.linspace(self.bandwith[0], self.bandwith[1], self.n)
        self.nu_arr = np.ones((self.n, self.n, self.n))*self.nu[:,None,None]


    #================FUNCTIONS RELATED TO THE REAL SPACE OUTPUT================#

    def generate_distances(self):
        '''
        The grid is centered at r_0, and the grid spacing is set by the real space
        length of the box, and the number of pixels along one axis. 

        The comoving distance of each slice is then determined by the resolution. 
        '''
        self.delta_r = self.L/self.n #realspace resolution
        self.comoving_dist = (np.ones(self.n)*self.r_0 
                - Distance((self.n/2 - np.arange(0,self.n))*self.delta_r, 
                unit = u.Mpc, allow_negative=True))


    def distances_to_freq(self):
        '''
        converts the comoving ditances of the grid back to frequencies.
        '''
        z = []
        for dist in self.comoving_dist:
            z.append(z_at_value(cosmo.comoving_distance, dist))
        self.redshifts = np.array(z) 

        self.frequencies = self.NU_21CM/(self.redshifts + 1)
        self.frequency_grid = np.ones((self.n, self.n, self.n))*self.frequencies[:,None,None]
        
    def freq_to_z(self, nu):
        '''
        Given a frequency nu, finds the redshift that would correspond to a 
        redshifted 21-cm photon detected with that frequency today. 
        '''
        return self.NU_21CM/nu -1

    def z_to_dist(self, z):
        ''' 
        Given a redshift, finds the correspoding comoving distance
        '''
        return cosmo.comoving_distance(z)


    #================FUNCTIONS RELATED TO THE FOREGROUND MODEL================#

   
    def generate_foregrounds(self, arr):
        '''generates total temperature per pixel due to all foregrounds'''

        #total temperature contributions to each pixel from all the foregrounds
        self.foreground_temps = np.zeros((self.n, self.n, self.n))
        self.amplitudes = []

        #one loop generates temperature due to one foreground in each pixel
        for i in range(self.n_f):
            self.foreground_temps += self.add_foreground_temp(arr)
        
        self.foreground_temps *= (1/self.n_f)

        self.amplitudes = np.array(self.amplitudes) #amplitudes of each foreground
   
    def add_foreground_temp(self, freq_array):
        '''for one foreground in each pixel, generates its spectrum 
        T(nu) = L*nu**(-alpha) where L is generated in generate_amplitudes() 
        and alpha is drawn from a normal distribution.
        
        Parameters
        ----------
        freq_array: 3D numpy array
            3D numpy array where one 2D slice has the frequencies at that index.
            Ie if nu[i] = nu_0: 
            freq_array[i] = [[nu_0,nu_0,....,nu_0],
                              ...
                             [nu_0,nu_0,....,nu_0]]
                             
        Returns
        -------
        temps: 3D numpy array
            a 3D numpy array where each pixel has temperature of one foreground 
            at a specific frequency
        '''
    
        self.generate_amplitudes()
        alpha = np.random.normal(self.f_alpha[0],self.f_alpha[1], (self.n,self.n))
        lnT = np.log(self.amps[None:,:]) - alpha*np.log(freq_array/150)
        temps = np.exp(lnT)
        return temps 

    def generate_amplitudes(self):
        '''The luminosity of each foreground, ie L in the foreground's spectrum
        T(nu) = L*nu**(-f_alpha) is drawn from a power law distribution 
        P(x) = N*x**(-l_alpha) for x > tmin, where N is a normalization factor.
        '''

        a = self.l_alpha -1 #since pareto is defined for x**(-1-a) not x**a
        
        m = self.mean_150Mhz*(a-1)/a #minimum temperature set by mean
        self.amps = (np.random.pareto(a, (self.n,self.n)) +1)*m 
        self.amplitudes.append(self.amps)



def main():

    f = Foregrounds()
    print(f.real_space())
    # print(f.r_0.value)

    test =Distance(20, unit = u.Mpc)
    # print(f.r_0)
    # print(f.delta_r)
    # print(f.comoving_dist)
    cdist = f.comoving_dist
    # z = []
    # for dist in cdist:
    #     z.append(z_at_value(cosmo.comoving_distance, dist))
    # # print(z_at_value(cosmo.comoving_distance, f.comoving_dist))
    # print(f.frequency_grid[100])
if __name__=="__main__":
    main()
