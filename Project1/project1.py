import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
import sys

class project1:
    def __init__(self, N, Omega_m0, Omega_lambda0, h):
        self.N = N
        self.Omega_m0 = Omega_m0
        self.Omega_lambda0 = Omega_lambda0
        self.h = h
        self.w = 1 -  Omega_m0 - Omega_lambda0
        self.c = 299792458    

    def integrate(self, x):
        """ Integral that solves for t*H_0"""
        return np.sqrt(x)/np.sqrt(self.Omega_m0 + self.w*x + self.Omega_lambda0*x**3)

    def S(self, r, w):
        """ The function Sk(r) which returns different depending on the value of r (w) """
        if w > 1:
            return np.sin(r)
        else:
            return np.sinh(r)


    def integrate_z(self, z, Omega_m0, Omega_lambda0, w):
        """ Integral that solves for t*H_0 """
        p = 1.0 + z
        return 1.0/np.sqrt(Omega_m0*p**3 + w*p**2 + Omega_lambda0)

    def integrate_time_z(self, z, Omega_m0, Omega_lambda0, w):
        """Integral that solves for t(z)"""
        p = 1.0 + z
        return 1.0/(p*np.sqrt(Omega_m0*p**3 + w*p**2 + Omega_lambda0))

    def find_age(self):
        """Method which returns t_0 and the value of H_0"""
        year = 60*60*24*365
        Mpc = 3.08567758e22
        km = 1e3    
        H_0 = 100*km*self.h/Mpc #[s-1]
        t_0 = quad(self.integrate, 0, 1)[0]/H_0 #[s-1]
        print "The age of the universe with current parameters is %.2f billion years" %(t_0/year/1e9)
        return H_0, t_0

    def m0_lambda0_plot(self):
        """Least square model to locate Omega-values within 5% range of empirical values """
        z = np.linspace(0, 10, self.N)              # liste of z-values
        M = 100                                      # nr. values in each list
        Omega_m0_list = np.linspace(0, 2, M)        
        Omega_lambda0_list = np.linspace(0, 2, M)
        tol = 1e-4                                  # tolerance for z-values
        tol_K = 0.02
        omegas_lambda0 = []
        omegas_m0 = []
        omegas_lambda0_K1 = []
        omegas_m0_K1 = []
        omegas_lambda0_K2 = []
        omegas_m0_K2 = []
        for i in range(M):
            Omega_m0 = Omega_m0_list[i]
            for j in range(M):
                dL = np.zeros(self.N)
                Omega_lambda0 = Omega_lambda0_list[j]
                w = 1 - Omega_lambda0 - Omega_m0
                for k in range(self.N):
                    self.test(Omega_m0, Omega_lambda0)  # checking the initial a/a_0 test
                    if exit == True:
                        break
                    if 0.0 < w <= 0.02:
                        omegas_lambda0_K1.append(Omega_lambda0)    # Saving value
                        omegas_m0_K1.append(Omega_m0)              # Saving value
                    elif -0.02 <= w < 0.00:
                        omegas_lambda0_K2.append(Omega_lambda0)    # Saving value
                        omegas_m0_K2.append(Omega_m0)              # Saving value
                    if z[k] > 0.84:                     # not interested in these values
                        break
                    elif w == 0:                        # k = 0
                        dL[k] = (1 + z[k])*quad(self.integrate_z, 0, z[k], args = (Omega_m0, Omega_lambda0, w))[0]
                    else:                               # using the S-function
                        dL[k]= (1 + z[k])/np.sqrt(abs(w))*self.S(np.sqrt(abs(w))*quad(self.integrate_z, 0, z[k], args = (Omega_m0, Omega_lambda0, w))[0], w)
                                
                    if abs(z[k] - 0.83) < tol:                      # Checking if we're close enough in the z-range
                        if 1.16*0.95 <= dL[k] <= 1.16*1.05:         # Inside the 5% margin
                            omegas_lambda0.append(Omega_lambda0)    # Saving value
                            omegas_m0.append(Omega_m0)              # Saving value
                            break
                        else:
                            pass        

        """
        M = 100                                      # nr. values in each list
        Omega_m0_list = np.linspace(0, 2, M)        
        Omega_lambda0_list = np.linspace(0, 2, M)
        omegas_lambda0_K1 = []
        omegas_m0_K1 = []
        omegas_lambda0_K2 = []
        omegas_m0_K2 = []

        for i in range(M):
            Omega_m0 = Omega_m0_list[i]
            for j in range(M):
                dL = np.zeros(self.N)
                Omega_lambda0 = Omega_lambda0_list[j]
                w = 1 - Omega_lambda0 - Omega_m0
                for k in range(self.N):
                    if 0.0 < w <= 0.02:
                        omegas_lambda0_K1.append(Omega_lambda0)    # Saving value
                        omegas_m0_K1.append(Omega_m0)              # Saving value
                        break
                    elif -0.02 <= w < 0.00:
                        omegas_lambda0_K2.append(Omega_lambda0)    # Saving value
                        omegas_m0_K2.append(Omega_m0)              # Saving value
                        break
        """
        plt.plot(omegas_m0, omegas_lambda0, "o")
        plt.ylabel("$\\Omega_{\\lambda_0}$")
        plt.xlabel("$\\Omega_{m_0}$")
        #plt.show()
        plt.show()

        plt.plot(omegas_m0_K1, omegas_lambda0_K1, "r-")
        plt.plot(omegas_m0_K2, omegas_lambda0_K2, "r-")
        plt.legend(["$\\Omega_{k0} \\pm 0.02$"])
        plt.plot(omegas_m0, omegas_lambda0, "o")
        plt.ylabel("$\\Omega_{\\lambda_0}$")
        plt.xlabel("$\\Omega_{m_0}$")
        #plt.show()
        plt.show()

    def test(self, Omega_m0, Omega_lambda0):
        x_test = np.linspace(0, 1, 100) # a_0/a test
        test = self.Omega_m0*x_test[1:]**(-3) + (1 - self.Omega_m0 - self.Omega_lambda0)*x_test[1:]**(-2) + self.Omega_lambda0
        for value in test:
            if value <= 0:
                print "The right hand side of the first Friedmann equation is not positive. Chose different parameters"    
                exit = True
            else:
                exit = False
    
    def run(self):        
        exit = False
        self.test(self.Omega_m0, self.Omega_lambda0)
        if exit == True:
            sys.exit()
        H_0, t_0 = self.find_age()
        x = np.zeros(self.N) # [a(t)/a_0]
        x[-1] = 1.0 #a(t_0)/a_0
        t = np.linspace(0, 1, self.N) #[in units of t/t_0]
        dt = t[1] - t[0]

        # tracing x backwards from x = 1 (a(t_o)/a_0) to x = a(t=0)/a_0
        for i in range(self.N-1, 0, -1):
            x[i-1] = x[i] - H_0*t_0*dt*np.sqrt(self.Omega_m0/x[i] + self.w + self.Omega_lambda0*x[i]**2)    
        plt.plot(t, x)
        
        if ((self.Omega_m0 == 1.0) and (self.Omega_lambda0 == 0.0)):
            plt.title("Einstein de Sitter model")

        else:
            plt.title("$Omega_0$ = %.1f and $Omega_{\\Lambda_0} = %.1f$" %(self.Omega_m0, self.Omega_lambda0))    
        
        plt.xlabel("$t/t_0$")
        plt.ylabel("$a/a_0$")
        plt.show()

        z = np.linspace(0, 10, self.N)
        time_z = np.zeros(self.N)
        for i in range(self.N):
           time_z[i] = quad(self.integrate_time_z, 0, z[i], args = (self.Omega_m0, self.Omega_lambda0, self.w))[0]/H_0

        plt.plot(z, time_z/t_0)
        plt.ylabel("$t/t_0$")
        plt.xlabel("z")
        plt.show()

        dL = np.zeros(self.N)
        for i in range(self.N):
            if self.w == 0:
                dL[i] = (1 + z[i])*quad(self.integrate_z, 0, z[i], args = (self.Omega_m0, self.Omega_lambda0, self.w))[0]
            else:
                dL[i]= (1 + z[i])/np.sqrt(abs(self.w))*self.S(np.sqrt(abs(self.w))*quad(self.integrate_z, 0, z[i], args = (self.Omega_m0, self.Omega_lambda0, self.w))[0], self.w)
        plt.plot(z, dL)
        plt.xlabel("Z")
        plt.ylabel("dL*c/H_0")
        plt.show()
        self.m0_lambda0_plot()

A = project1(int(1e4), 1.0, 0.0, 0.72)
A.run()
