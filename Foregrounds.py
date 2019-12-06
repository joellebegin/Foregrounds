import numpy as np 

class Foregrounds:
    ''' generates 3d box with extragalactic foregrounds. Each pixel has a
    specified number of foregrounds, where every foreground has a different 
    spectrum and luminosity as specified in generate_amplitudes() and
    add_foreground_temp().

    -bandwith: range of frequencies that will be considered in the parallel 
    direction. Array/tuple with [nu_min, nu_max]
    -n: dimensions of box
    -foregrounds_per_pixel: number of foregrounds in one pixel'''


    def generate_amplitudes(self):
        '''The luminosity of each foreground in the first frequency slice is 
        drawn from a power law distribution P(x) = a*(m**a)*x**(-1-a) where m is
        the minimum temperature to be considered.
        '''
        a = 3
        m = 2
        self.amps = (np.random.pareto(a, (self.n,self.n)) +1)*m 
        
        self.amplitudes += self.amps

    def add_foreground_temp(self):
        '''for each foreground, generates its spectrum T(nu) = A*nu**(-alpha)
        where A is generated in generate_amplitudes() and alpha is drawn from 
        a normal distribution.'''

        self.generate_amplitudes()
        alpha = np.random.normal(2.5,0.5, (self.n,self.n))
        lnT = np.log(self.amps[None:,:]) - alpha*np.log(self.nu_arr)
        self.T += np.exp(lnT)

    def generate_foregrounds(self):
        '''generates total temperature per pixel due to all foregrounds'''
        
        #total temperature contributions to each pixel from all the foregrounds
        self.T = np.zeros((self.n-1, self.n, self.n))
        self.amplitudes = np.zeros((self.n,self.n))
        
        for i in range(self.n_f):
            #one loop generates temperature due to one foreground in each pixel
            self.add_foreground_temp()
        
        #first frequency slice is sum of luminosities generated in generate_amplitudes()
        self.temp = np.vstack((self.amplitudes[None], self.T))


    def generate_frequency_array(self):

        self.nu = np.linspace(self.bandwith[0], self.bandwith[1], self.n -1)
        self.nu_arr = np.ones((self.n-1, self.n, self.n))*self.nu[:,None,None]
        
    def __init__(self, bandwith, n, foregrounds_per_pix):
        
        self.n = n 
        self.n_f = foregrounds_per_pix
        self.bandwith = bandwith

        
        self.generate_frequency_array()
        self.generate_foregrounds()

def main():

    f = Foregrounds([90,100], 4, 1)
    print(f.temp)

if __name__=="__main__":
    main()