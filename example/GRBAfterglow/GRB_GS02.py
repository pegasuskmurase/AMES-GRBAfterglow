import numpy as np

class GS02:
    """
    GRB afterglow lightcurve and spectrum from Granot & Sari 2002.
    """
    def __init__(self, E, n, epsilon_e, epsilon_B, p, dl, z):
        super().__init__()
        self.E = E
        self.n = n
        self.epsilon_e = epsilon_e
        self.epsilon_e_bar = epsilon_e * (p-2) / (p-1)
        self.epsilon_B = epsilon_B
        self.p = p
        self.dl = dl
        self.z = z

    def gen_lightcurve_data(self, nu_l, tdays_l):
        beta1_l = [0, 2, 1./3, (1-self.p)/2, 2, 5./2]
        beta2_l = [0, 1./3, (1-self.p)/2, -self.p/2, 5./2, (1-self.p)/2]
        s_l = [0, 1.64, 1.84-0.4*self.p, 1.15 - 0.06*self.p, 3.44*self.p-1.41, 1.47 - 0.21*self.p]

        flux = []
        for nu, tdays in zip(nu_l, tdays_l):
            b = 1
            f1 = self.f_density(nu, self.calc_nu_b(b, tdays), self.calc_f_nu_b(b, tdays), beta1_l[b], beta2_l[b], s_l[b])
            b = 2
            f2 = self.f_density_wave(nu, self.calc_nu_b(b, tdays), beta1_l[b], beta2_l[b], s_l[b])
            b = 3
            f3 = self.f_density_wave(nu, self.calc_nu_b(b, tdays), beta1_l[b], beta2_l[b], s_l[b])
            b = 4
            f4 = self.f_density4(nu, self.calc_nu_b(b, tdays), self.calc_f_nu_b(b, tdays), beta1_l[b], beta2_l[b], s_l[b])
            b = 5
            f5 = self.f_density_wave(nu, self.calc_nu_b(b, tdays), beta1_l[b], beta2_l[b], s_l[b])

            nu1 = self.calc_nu_b(1, tdays)
            nu2 = self.calc_nu_b(2, tdays)

            if nu1 < nu2:
                flux.append(f1*f2*f3*1e3)
            elif nu1 >= nu2:
                flux.append(f4*f5*f3*1e3)
            else:
                flux.append(0)
                print('out of range !!!')
        return np.array(flux)

    def gen_lightcurve(self, nu, tdays_l):
        beta1_l = [0, 2, 1./3, (1-self.p)/2, 2, 5./2]
        beta2_l = [0, 1./3, (1-self.p)/2, -self.p/2, 5./2, (1-self.p)/2]
        s_l = [0, 1.64, 1.84-0.4*self.p, 1.15 - 0.06*self.p, 3.44*self.p-1.41, 1.47 - 0.21*self.p]

        flux = []
        for tdays in tdays_l:
            b = 1
            f1 = self.f_density(nu, self.calc_nu_b(b, tdays), self.calc_f_nu_b(b, tdays), beta1_l[b], beta2_l[b], s_l[b])
            b = 2
            f2 = self.f_density_wave(nu, self.calc_nu_b(b, tdays), beta1_l[b], beta2_l[b], s_l[b])
            b = 3
            f3 = self.f_density_wave(nu, self.calc_nu_b(b, tdays), beta1_l[b], beta2_l[b], s_l[b])
            b = 4
            f4 = self.f_density4(nu, self.calc_nu_b(b, tdays), self.calc_f_nu_b(b, tdays), beta1_l[b], beta2_l[b], s_l[b])
            b = 5
            f5 = self.f_density_wave(nu, self.calc_nu_b(b, tdays), beta1_l[b], beta2_l[b], s_l[b])
            f1 

            nu1 = self.calc_nu_b(1, tdays)
            nu2 = self.calc_nu_b(2, tdays)

            if nu1 < nu2:
                flux.append(f1*f2*f3)
            elif nu1 >= nu2:
                flux.append(f4*f5*f3)
            else:
                print('out of range !!!')
        return np.array(flux)

    def gen_spectrum(self, nu, tdays):
        beta1_l = [0, 2, 1./3, (1-self.p)/2, 2, 5./2]
        beta2_l = [0, 1./3, (1-self.p)/2, -self.p/2, 5./2, (1-self.p)/2]
        s_l = [0, 1.64, 1.84-0.4*self.p, 1.15 - 0.06*self.p, 3.44*self.p-1.41, 1.47 - 0.21*self.p]

        b = 1
        f1 = self.f_density(nu, self.calc_nu_b(b, tdays), self.calc_f_nu_b(b, tdays), beta1_l[b], beta2_l[b], s_l[b])
        f1_norm = self.f_density(1e14, self.calc_nu_b(b, tdays), self.calc_f_nu_b(b, tdays), beta1_l[b], beta2_l[b], s_l[b])
        b = 2
        f2 = self.f_density_wave(nu, self.calc_nu_b(b, tdays), beta1_l[b], beta2_l[b], s_l[b])
        f2_norm = self.f_density_wave(1e14, self.calc_nu_b(b, tdays), beta1_l[b], beta2_l[b], s_l[b])
        b = 3
        f3 = self.f_density_wave(nu, self.calc_nu_b(b, tdays), beta1_l[b], beta2_l[b], s_l[b])
        f3_norm = self.f_density_wave(1e14, self.calc_nu_b(b, tdays), beta1_l[b], beta2_l[b], s_l[b])
        b = 4
        f4 = self.f_density4(nu, self.calc_nu_b(b, tdays), self.calc_f_nu_b(b, tdays), beta1_l[b], beta2_l[b], s_l[b])
        f4_norm = self.f_density(1e14, self.calc_nu_b(b, tdays), self.calc_f_nu_b(b, tdays), beta1_l[b], beta2_l[b], s_l[b])
        b = 5
        f5 = self.f_density_wave(nu, self.calc_nu_b(b, tdays), beta1_l[b], beta2_l[b], s_l[b])
        f5_norm = self.f_density_wave(1e14, self.calc_nu_b(b, tdays), beta1_l[b], beta2_l[b], s_l[b])


        nu1 = self.calc_nu_b(1, tdays)
        nu2 = self.calc_nu_b(2, tdays)
        nu3 = self.calc_nu_b(3, tdays)

        Hz2eV = 4.1e-15
        print('GS02 eps_a = ', nu1*Hz2eV, 'eps_m = ', nu2*Hz2eV, 'eps_c = ', nu3*Hz2eV)

        if nu1 < nu2:
            #print('norm = ', self.normalization('B', tdays))
            #norm = self.normalization('B', tdays)/f1_norm/f2_norm/f3_norm
            #print(norm)
            return f1*f2*f3
        elif nu1 >= nu2:
            print(tdays)
            #norm = self.normalization('B', tdays)/f4_norm/f5_norm/f3_norm
            return f4*f5*f3
        else:
            print('out of range !!!')


    def normalization(self, PLS, tdays):
        nu_norm = 1e14
        if PLS == 'B':
            return 4.2 * (3*self.p+2) / (3*self.p-1) * 1e9 * (1+self.z)**(5./2) * self.epsilon_e_bar * self.n**(-1./2) * (self.E/1e52)**(1./2) * tdays**(1./2) * (self.dl/1e28)**(-2) * (nu_norm/1e14)**2
        else:
            print("out of range!!!")

    def f_density_wave(self, nu, nu_b, beta1, beta2, s):
        return (1 + (nu/nu_b)**(s*(beta1 - beta2)))**(-1/s)

    def f_density(self, nu, nu_b, f_nu_b, beta1, beta2, s):
        return f_nu_b * ((nu/nu_b)**(-s*beta1) + (nu/nu_b)**(-s*beta2))**(-1/s)

    def f_density4(self, nu, nu_b, f_nu_b, beta1, beta2, s):
        phi4 = nu / nu_b
        return f_nu_b * (phi4**2 * np.exp(-s*phi4**(2./3)) + phi4**(5./2))

    def calc_nu_b(self, b, tdays):
        if b == 1:
            return 1.24 * (self.p-1)**(3./5) / (3*self.p+2.)**(3./5) * 1e9 * (1+self.z)**(-1) * self.epsilon_e_bar**(-1) * self.epsilon_B**(1/5) * self.n**(3./5) * (self.E/1e52)**(1./5)
        elif b == 2:
            return 3.73 * (self.p - 0.67) * 1e15 * (1 + self.z)**(1./2) * (self.E/1e52)**(1./2) * self.epsilon_e_bar**(2) * self.epsilon_B**(1/2) * tdays**(-3./2)
        elif b == 3:
            return 6.37 * (self.p-0.46) * 1e13 * np.exp(-1.16*self.p) * (1 + self.z)**(-1./2) * self.epsilon_B**(-3./2) * self.n**(-1) * (self.E/1e52)**(-1./2) * tdays**(-1./2)
        elif b == 4:
            return 5.04 * (self.p-1.22) * 1e16 * (1 + self.z)**(1./2) * self.epsilon_e_bar**(2) * self.epsilon_B**(1./2) * (self.E/1e52)**(1./2) * tdays**(-3./2)
        elif b == 5:
            return 3.59 * (4.03 - self.p) * 1e9 * np.exp(2.34*self.p) * (self.epsilon_e_bar**(4*(self.p-1)) / (1 + self.z)**(6-self.p) * self.epsilon_B**(self.p + 2) * self.n**4 * (self.E/1e52)**(self.p+2) / tdays**(self.p*3+2))**(1./2/(self.p+4))
        else:
            print('out of range!!!')

    def calc_f_nu_b(self, b, tdays):
        if b == 1:
            return 0.647 * (self.p-1)**(6./5) / (3*self.p+2.)**(1./5) / (3*self.p-1) * (1+self.z)**(1./2) * self.epsilon_e_bar**(-1) * self.epsilon_B**(2./5) * self.n**(7./10) * (self.E/1e52)**(9./10) * tdays**(1./2) * (self.dl/1e28)**(-2) 
        elif b == 2:
            return 9.93 * (self.p + 0.14) * (1 + self.z) * (self.E/1e52) * self.epsilon_B**(1./2) * self.n**(1/2) * (self.dl/1e28)**(-2)
        elif b == 3:
            return 4.68 * np.exp(4.82*(self.p-2.5)) * 1e3 * (1 + self.z)**((self.p+1)/2) * self.epsilon_e_bar**(self.p-1) * self.epsilon_B**(self.p-1./2) * self.n**(self.p/2) * (self.E/1e52)**((self.p+1)/2) * tdays**((1-self.p)/2) * (self.dl/1e28)**(-2)
        elif b == 4:
            return 3.72 * (self.p-1.79) * 1e15 * (1 + self.z)**(7./2) * self.epsilon_e_bar**5 * self.epsilon_B * self.n**(-1./2) * (self.E/1e52)**(3./2) * tdays**(-5./2) * (self.dl/1e28)**(-2)
        elif b == 5:
            return 20.8 * (self.p-1.53) * np.exp(2.56*self.p) * (self.dl/1e28)**(-2)  * ((1 + self.z)**(7*self.p+3) / self.epsilon_e_bar**(10.*(1-self.p)) * self.epsilon_B**(2*self.p+3) * (self.E/1e52)**(3*self.p+7) / tdays**(5*(self.p-1)))**(1./2/(self.p+4)) 
        else:
            print('out of range!!!')

