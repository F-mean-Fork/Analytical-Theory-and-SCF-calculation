from ascf import FlatBrush, TwoBrush
import matplotlib.pyplot as plt
import numpy as np
import warnings
warnings.filterwarnings('ignore')

if __name__ == "__main__":  
    params = [[112, 8, 1], [256, 24, 8], [1000, 0, 10]]
    colors = {3.0: 'red', 2.0: 'blue', 1.0: 'green'}
    for p in params:
        Nb, n, m = p 
        N = 1000
        eta = (n/m+1)**0.5
        chi = 0.0
        sigma = 0.01
        tb = TwoBrush(N, sigma, chi, eta)
        X = np.loadtxt(f"Nb{Nb}n{n}m{m}sigma{sigma}chi{chi}.txt")
        D = X[0]
        L = X[1][1::2]
        G = X[2][1::2]
        Phi = X[3][1::2]
        plt.scatter(Phi**2*L, G, color=colors[eta], s=12, marker='o', label=f"$\eta={eta}$")

    plt.ylabel("$\Omega$")
    plt.xlabel("$\phi_m^2$")
    plt.legend()
    plt.show()

  
    # Gamma on D
    params = [[112, 8, 1], [256, 24, 8], [1000, 0, 10]]
    colors = {3.0: 'red', 2.0: 'blue', 1.0: 'green'}
    for p in params:
        Nb, n, m = p 
        N = 1000
        eta = (n/m+1)**0.5
        chi = 0.0
        sigma = 0.01
        tb = TwoBrush(N, sigma, chi, eta)
        X = np.loadtxt(f"Nb{Nb}n{n}m{m}sigma{sigma}chi{chi}.txt")
        D = X[0]
        L = X[1]
        G = X[2][0::2]
        plt.scatter(D[0::2], G, color=colors[eta], s=12, marker='o', label=f"$\eta={eta}$")
        plt.plot(2*tb.d, [tb.overlap_integral(0.2, d) for d in tb.d], '--', color=colors[eta])
    plt.loglog()
    plt.ylim(1e-3, 5)
    plt.xlabel('$D$', fontsize=14)
    plt.ylabel('$\Gamma(D)$', fontsize=14)
    plt.title("$N = 1000, \sigma = 0.01, \chi = 0.0$")
    plt.legend(fontsize=14)
    plt.savefig('Gamma_on_D.png', bbox_inches='tight')
    plt.close()
    
    # Pi on D
    params = [[112, 8, 1], [256, 24, 8], [1000, 0, 10]]
    colors = {3.0: 'red', 2.0: 'blue', 1.0: 'green'}
    for p in params:
        Nb, n, m = p 
        N = 1000
        eta = (n/m+1)**0.5
        chi = 0.0
        sigma = 0.01
        tb = TwoBrush(N, sigma, chi, eta)
        
        X = np.loadtxt(f"Nb{Nb}n{n}m{m}sigma{sigma}chi{chi}.txt")
        D = X[0]
        L = X[1]
        F = X[4]
        Pi = -np.diff(F)/np.diff(D)
        D = 0.5*(D[1:] + D[:-1])
        plt.scatter(D[0:-1:5], Pi[0:-1:5], color=colors[eta], s=12,marker='o', label=f"$\eta={eta}$")
        plt.plot(2*tb.d, [tb.pressure(d) for d in tb.d], '--', color=colors[eta])
    plt.loglog()
    plt.ylim(1e-3, 3)
    plt.xlabel('$D$', fontsize=14)
    plt.ylabel('$\Pi(D)$', fontsize=14)
    plt.title("$N = 1000, \sigma = 0.01, \chi = 0.0$")
    plt.legend(fontsize=14)
    plt.savefig('Pi_on_D.png', bbox_inches='tight')
    plt.close()
        
        
    # f on Pi
    k = (np.pi - 2) / 4
    params = [[112, 8, 1], [256, 24, 8], [1000, 0, 10]]
    colors = {3.0: 'red', 2.0: 'blue', 1.0: 'green'}
    for p in params:
        Nb, n, m = p 
        N = 1000
        eta = (n/m+1)**0.5
        chi = 0.0
        sigma = 0.01
        tb = TwoBrush(N, sigma, chi, eta)
        
        X = np.loadtxt(f"Nb{Nb}n{n}m{m}sigma{sigma}chi{chi}.txt")
        D = X[0]
        L = X[1]
        G = X[2]*k           #OmegaL
        F = X[4]
        Pi = -np.diff(F)/np.diff(D)
        D = 0.5*(D[1:] + D[:-1])
        G = 0.5*(G[1:] + G[:-1])
        plt.scatter(Pi[0:-1:5], G[0:-1:5], color=colors[eta], s=12,marker='o', label=f"$\eta={eta}$")
        plt.plot([tb.pressure(d) for d in tb.d], [k*tb.overlap_integral(0.2, d) for d in tb.d], '--', color=colors[eta])
    plt.loglog()
    plt.ylim(1e-3, 3)
    plt.xlim(1e-3, 3)
    plt.xlabel('$\Pi$', fontsize=14)
    plt.ylabel('$f(\Pi)/(V \zeta)$', fontsize=14)
    plt.title("$N = 1000, \sigma = 0.01, \chi = 0.0$")
    plt.legend(fontsize=14)
    plt.savefig('f_on_Pi.png', bbox_inches='tight')
    plt.close()
    
    # # mu on Pi
    k = (np.pi - 2) / 4
    params = [[112, 8, 1], [256, 24, 8], [1000, 0, 10]]
    colors = {3.0: 'red', 2.0: 'blue', 1.0: 'green'}
    for p in params:
        Nb, n, m = p 
        N = 1000
        eta = (n/m+1)**0.5
        chi = 0.0
        sigma = 0.01
        tb = TwoBrush(N, sigma, chi, eta)
        
        X = np.loadtxt(f"Nb{Nb}n{n}m{m}sigma{sigma}chi{chi}.txt")
        D = X[0]
        L = X[1]
        G = X[2]*k           #OmegaL
        F = X[4]
        Pi = -np.diff(F)/np.diff(D)
        D = 0.5*(D[1:] + D[:-1])
        G = 0.5*(G[1:] + G[:-1])
        plt.scatter(Pi[0:-1:5], G[0:-1:5]/Pi[0:-1:5], color=colors[eta], s=12, marker='o', label=f"$\eta={eta}$")
        plt.plot([tb.pressure(d) for d in tb.d], [k*tb.overlap_integral(0.2, d)/tb.pressure(d) for d in tb.d], '--', color=colors[eta])
    plt.loglog()
    plt.ylim(0.0, 3.0)
    plt.xlim(1e-3, 3)
    plt.xlabel('$\Pi$', fontsize=14)
    plt.ylabel('$\mu/(V \zeta)$', fontsize=14)
    plt.title("$N = 1000, \sigma = 0.01, \chi = 0.0$")
    plt.legend(fontsize=14)
    plt.savefig('mu_on_Pi.png', bbox_inches='tight')
    plt.close()
