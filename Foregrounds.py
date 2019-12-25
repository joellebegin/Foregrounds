import numpy as np 
from astropy import units as u 
from astropy.coordinates import Distance

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

    def generate_frequency_array(self):
        '''generates the range of frequencies considered. Max and min frequencies
        are given by the specified bandwith, and the spacing is given by the
        size of the box, n.'''

        self.nu = np.linspace(self.bandwith[0], self.bandwith[1], self.n)
        self.nu_arr = np.ones((self.n, self.n, self.n))*self.nu[:,None,None]
    

    def generate_foregrounds(self):
        '''generates total temperature per pixel due to all foregrounds, where 
        the main axis is the frequency axis. i.e. frequency_space_temps[i] gives
        the temperatures at a frequency of self.nu[i]'''
        
        #total temperature contributions to each pixel from all the foregrounds
        self.frequency_space_temps = np.zeros((self.n, self.n, self.n))
        
        #one loop generates temperature due to one foreground in each pixel
        for i in range(self.n_f):
            self.add_foreground_temp()
        
    def add_foreground_temp(self):
        '''for each foreground, generates its spectrum T(nu) = L*nu**(-alpha)
        where L is generated in generate_amplitudes() and alpha is drawn from 
        a normal distribution.'''

        self.generate_amplitudes()
        alpha = np.random.normal(self.f_alpha[0],self.f_alpha[1], (self.n,self.n))
        lnT = np.log(self.amps[None:,:]) - alpha*np.log(self.nu_arr)
        self.frequency_space_temps += np.exp(lnT) 

    def generate_amplitudes(self):
        '''The luminosity of each foreground, ie L in the foreground's spectrum
        T(nu) = L*nu**(-f_alpha) is drawn from a power law distribution 
        P(x) = N*x**(-l_alpha) for x > tmin, where N is a normalization factor.
        '''

        a = self.l_alpha -1 #since pareto is defined for x**(-1-a) not x**a
        m = self.tmin
        self.amps = (np.random.pareto(a, (self.n,self.n)) +1)*m 
        
    def realspace(self):
        '''converts the frequency space axis into real space distance'''
        
        self.freq_to_z() #converts nu --> z
        self.z_to_dist() #converts z --> r

    def freq_to_z(self):
        '''for each frequency in self.nu, converts to the redshift corresponding 
        to a 21-cm photon'''

        self.z_arr = self.NU_21CM/self.nu -1

    def z_to_dist(self):
        '''for each redshift in self.z_arr, converts to the distance that that
        corresponds to'''

        self.dist_arr = Distance(unit = u.pc, z = self.z_arr)


def main():

    f = Foregrounds()
    freqs = f.frequency_space()

if __name__=="__main__":
    main()