import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
#Python Module to develop digital signal processing. It analyses an array of data, plotting it in time and frequency domain. It is able to plot
#the result of a digital filtering and compare it with the original one.

class DSP:
    def __init__(self, data, fs, filt = 'No', num = np.array([1]), den = np.array([1, 0])):
        self.data = np.array(data)
        self.filt = filt
        self.num = num
        self.den = den
        self.data_l = list(data)
        self.data_f = np.array([])
        self.fs = fs
        self.n = len(data)
        self.F = fs/(self.n)
        self.YFFT = []
        self.YF_FFT = []
        self.x = []
        self.XFFT = np.array(range(0,round(self.n/2)))
        self.XFFT = self.F*self.XFFT
        for i in list(self.XFFT):
            self.x.append(i)
        
    #It uses the FFT algorithm to transform discrete time domain data in discrete frequency domain
    def calculate(self,):
        YFFT2 = abs(np.fft.fft(self.data))/len(self.data)
        n = 0
        for i in list(YFFT2):
            if (n == 0):
                self.YFFT.append(i)
            elif (n < round(self.n/2)):
                self.YFFT.append(2*i)
            else:
                break
            n = n + 1
        #If a filter has been defined, it calculates how the data array would have been using this filter
        if (self.filt == 'Yes'):
            y = signal.filtfilt(self.num, self.den, self.data_l)
            self.data_f = y
            YFFT3 = abs(np.fft.fft(y))/len(y)
            n = 0
            for i in list(YFFT3):
                if (n == 0):
                    self.YF_FFT.append(i)
                elif (n < round(self.n/2)):
                    self.YF_FFT.append(2*i)
                else:
                    break
                n = n + 1
    
    #Once the numerical calculations have been done, they can be plotted to see the results
    def plot(self,):
        t = list(np.linspace(0,1/self.F,self.n))
        fig, ax = plt.subplots(2)
        fig.suptitle('Frequencial Analysis of the Signal')
        ax[0].plot(t, self.data_l,'b', label = 'Original')
        if (self.filt == 'Yes'):
            ax[0].plot(t, self.data_f,'r', label = 'Filtered')
        ax[0].set(xlabel='Seconds', ylabel='Value')
        ax[0].grid()
        ax[0].legend(loc='best')
        ax[1].plot(self.x, self.YFFT,'b', label = 'Original')
        if (self.filt == 'Yes'):
            ax[1].plot(self.x, self.YF_FFT,'r', label = 'Filtered')
        ax[1].grid()
        ax[1].axis([min(self.x)-0.01*max(self.x), max(self.x), 0-0.01*max(self.YFFT), 1.05*max(self.YFFT)])
        ax[1].set(xlabel='Frequency (Hz)', ylabel='Amplitude')
        ax[1].legend(loc='best')
        plt.show()
        