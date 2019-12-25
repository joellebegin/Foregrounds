import numpy as np 
from astropy.cosmology import WMAP9 as cosmo

class Foregrounds:
    ''' generates 3d box with extragalactic foregrounds. Each pixel has a
    specified number of foregrounds, where every foreground has a different 
    spectrum T(nu) = A*v**(-alpha), where A is as specified in generate_amplitudes().

    -n: dimensions of box
    -foregrounds_per_pixel: number of foregrounds in one pixel
    -foreground_alpha: specifies the gaussian distribution from which the index 
    of the foreground is drawn. Array/tuple with [mean, std]
    -min_temp: mininum temperature to be considered for the luminosity of foregrounds
    -luminosty_alpha: index describing the power law distribution from which foreground
    luminosities will be drawn.'''

    def __init__(self, n=200, foregrounds_per_pix=5, foreground_alpha = [2.5,0.5], 
                min_temp = 3., luminosity_alpha = 4.):
        
        self.n = n 
        self.n_f = foregrounds_per_pix
        self.f_alpha = foreground_alpha
        self.tmin = min_temp
        self.l_alpha = luminosity_alpha
        self.NU_21CM = 1420 #MHz


    def frequency_space(self, bandwith = [100,200]):
        '''-bandwith: range of frequencies that will be considered. 
        Array/tuple with [nu_min, nu_max]'''

        self.bandwith = bandwith

        #foregrounds in frequency space
        self.generate_frequency_array()
        self.generate_foregrounds()

        return(self.frequency_space_temps)

    def real_space(self, central_freq = 150, L = 300):
        
        self.r_0 = self.z_to_dist(self.freq_to_z(central_freq))


    def generate_frequency_array(self): #fspace
        '''generates the range of frequencies considered. Max and min frequencies
        are given by the specified bandwith, and the spacing is given by the
        size of the box, n.'''

        self.nu = np.linspace(self.bandwith[0], self.bandwith[1], self.n)
        self.nu_arr = np.ones((self.n, self.n, self.n))*self.nu[:,None,None]
    

    def generate_foregrounds(self): #fspace
        '''generates total temperature per pixel due to all foregrounds, where 
        the main axis is the frequency axis. i.e. frequency_space_temps[i] gives
        the temperatures at a frequency of self.nu[i]'''
        
        #total temperature contributions to each pixel from all the foregrounds
        self.frequency_space_temps = np.zeros((self.n, self.n, self.n))
        
        #one loop generates temperature due to one foreground in each pixel
        for i in range(self.n_f):
            self.frequency_space_temps += self.add_foreground_temp(self.nu_arr)
        

    def add_foreground_temp(self, freq_array): #foregrounds
        '''for each foreground, generates its spectrum T(nu) = L*nu**(-alpha)
        where L is generated in generate_amplitudes() and alpha is drawn from 
        a normal distribution.'''

        self.generate_amplitudes()
        alpha = np.random.normal(self.f_alpha[0],self.f_alpha[1], (self.n,self.n))
        lnT = np.log(self.amps[None:,:]) - alpha*np.log(freq_array)
        return np.exp(lnT) 

    def generate_amplitudes(self): #foregrounds
        '''The luminosity of each foreground, ie L in the foreground's spectrum
        T(nu) = L*nu**(-f_alpha) is drawn from a power law distribution 
        P(x) = N*x**(-l_alpha) for x > tmin, where N is a normalization factor.
        '''

        a = self.l_alpha -1 #since pareto is defined for x**(-1-a) not x**a
        m = self.tmin
        self.amps = (np.random.pareto(a, (self.n,self.n)) +1)*m 

    def freq_to_z(self, nu):
    
        return self.NU_21CM/nu -1

    def z_to_dist(self, z):
       
        return cosmo.comoving_distance(z)


def main():

    f = Foregrounds()
    freqs = f.frequency_space()

if __name__=="__main__":
    main()