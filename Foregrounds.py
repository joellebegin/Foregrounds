import numpy as np 
from astropy.cosmology import WMAP9 as cosmo

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

    min_temp: float 
        mininum temperature to be considered for the luminosity of foregrounds. 
        Very important that it is float and not int, the pareto function does not
        like ints. 

    luminosty_alpha: float 
        index describing the power law distribution from which foreground
        luminosities will be drawn. Also very important that float and not int

    
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
                min_temp = 3., luminosity_alpha = 4.):
        
        self.n = n 
        self.n_f = foregrounds_per_pix
        self.f_alpha = foreground_alpha
        self.tmin = min_temp
        self.l_alpha = luminosity_alpha
        self.NU_21CM = 1420 #MHz


    #========================METHODS THAT OUTPUT BOXES========================#

    def frequency_space(self, bandwith = [100,200]):
        '''
        Returns a box that is uniformly spaced in frequency. The range of frequencies
        is specified by the bandwidth, and the spacing by the size of the box.
        Every slice contains the foreground temperature at a specific frequency. 
        That is, frequency_space_temps[i] is what the sky would look like at 
        a frequency of self.nu[i]. 

        Parameters
        ----------
        -bandwith: Array/tuple of ints or floats
            range of frequencies that will be considered. 
            Array/tuple with [nu_min, nu_max]
            
        Returns
        -------
        -frequency_space_temps: 3D numpy array
            3D box where every slice contains foreground temperatures at different
            frequencies
        -nu: 1D numpy array
            The frequencies corresponding to each slice. frequency_space_temps[i]
            is at frequency nu[i]
        '''

        self.bandwith = bandwith

        #foregrounds in frequency space
        self.generate_frequency_array()
        self.generate_foregrounds()

        return(self.frequency_space_temps)

    def real_space(self, central_freq = 150, L = 300):
        
        self.r_0 = self.z_to_dist(self.freq_to_z(central_freq))


    #=============FUNCTIONS RELATED TO THE FREQUENCY SPACE OUTPUT=============#

    def generate_frequency_array(self):
        '''generates the range of frequencies considered.'''

        self.nu = np.linspace(self.bandwith[0], self.bandwith[1], self.n)
        self.nu_arr = np.ones((self.n, self.n, self.n))*self.nu[:,None,None]
    

    def generate_foregrounds(self):
        '''generates total temperature per pixel due to all foregrounds'''

        #total temperature contributions to each pixel from all the foregrounds
        self.frequency_space_temps = np.zeros((self.n, self.n, self.n))
        
        #one loop generates temperature due to one foreground in each pixel
        for i in range(self.n_f):
            self.frequency_space_temps += self.add_foreground_temp(self.nu_arr)
        
    
    #================FUNCTIONS RELATED TO THE REAL SPACE OUTPUT================#

    def freq_to_z(self, nu):
    
        return self.NU_21CM/nu -1

    def z_to_dist(self, z):
       
        return cosmo.comoving_distance(z)


    #================FUNCTIONS RELATED TO THE FOREGROUND MODEL================#

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
        lnT = np.log(self.amps[None:,:]) - alpha*np.log(freq_array)
        temps = np.exp(lnT)
        return temps 

    def generate_amplitudes(self):
        '''The luminosity of each foreground, ie L in the foreground's spectrum
        T(nu) = L*nu**(-f_alpha) is drawn from a power law distribution 
        P(x) = N*x**(-l_alpha) for x > tmin, where N is a normalization factor.
        '''

        a = self.l_alpha -1 #since pareto is defined for x**(-1-a) not x**a
        m = self.tmin
        self.amps = (np.random.pareto(a, (self.n,self.n)) +1)*m 



def main():

    f = Foregrounds()
    freqs = f.frequency_space()

if __name__=="__main__":
    main()