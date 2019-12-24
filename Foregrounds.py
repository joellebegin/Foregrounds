import numpy as np 

class Foregrounds:
    ''' generates 3d box with extragalactic foregrounds. Each pixel has a
    specified number of foregrounds, where every foreground has a different 
    spectrum T(nu) = A*v**(-alpha), where A is as specified in generate_amplitudes().

    -bandwith: range of frequencies that will be considered. Array/tuple with [nu_min, nu_max]
    -n: dimensions of box
    -foregrounds_per_pixel: number of foregrounds in one pixel
    -foreground_alpha: specifies the gaussian distribution from which the index 
    of the foreground is drawn. Array/tuple with [mean, std]
    -min_temp: mininum temperature to be considered for the luminosity of foregrounds
    -luminosty_alpha: index describing the power law distribution from which foreground
    luminosities will be drawn.'''

    def __init__(self, bandwith, n, foregrounds_per_pix, foreground_alpha = [2.5,0.5],
                min_temp = 3., luminosity_alpha = 4.):
        
        self.n = n 
        self.n_f = foregrounds_per_pix
        self.bandwith = bandwith
        self.f_alpha = foreground_alpha
        self.tmin = min_temp
        self.l_alpha = luminosity_alpha

        self.generate_frequency_array()
        self.generate_foregrounds()

    def generate_frequency_array(self):
        '''generates the range of frequencies considered. Max and min frequencies
        are given by the specified bandwith, and the spacing is given by the
        size of the box, n.'''

        self.nu = np.linspace(self.bandwith[0], self.bandwith[1], self.n -1)
        self.nu_arr = np.ones((self.n-1, self.n, self.n))*self.nu[:,None,None]
    

    def generate_foregrounds(self):
        '''generates total temperature per pixel due to all foregrounds'''
        
        #total temperature contributions to each pixel from all the foregrounds
        self.T = np.zeros((self.n-1, self.n, self.n))
        self.amplitudes = np.zeros((self.n,self.n))
        
        #one loop generates temperature due to one foreground in each pixel
        for i in range(self.n_f):
            self.add_foreground_temp()
        
        #first frequency slice is sum of luminosities generated in generate_amplitudes()
        self.temp = np.vstack((self.amplitudes[None], self.T))

    def add_foreground_temp(self):
        '''for each foreground, generates its spectrum T(nu) = L*nu**(-alpha)
        where L is generated in generate_amplitudes() and alpha is drawn from 
        a normal distribution.'''

        self.generate_amplitudes()
        alpha = np.random.normal(self.f_alpha[0],self.f_alpha[1], (self.n,self.n))
        lnT = np.log(self.amps[None:,:]) - alpha*np.log(self.nu_arr)
        self.T += np.exp(lnT) 

    def generate_amplitudes(self):
        '''The luminosity of each foreground, ie L in the foreground's spectrum
        T(nu) = L*nu**(-f_alpha) is drawn from a power law distribution 
        P(x) = N*x**(l_alpha) for x > tmin, where N is a normalization factor.
        '''

        a = self.l_alpha -1 #since paretor is defined for x**(-1-a) not x**a
        m = self.tmin
        self.amps = (np.random.pareto(a, (self.n,self.n)) +1)*m 
        
        self.amplitudes += self.amps


def main():

    f = Foregrounds([90,100], 4, 1)
    print(f.nu)
    print(f.nu_arr)

if __name__=="__main__":
    main()