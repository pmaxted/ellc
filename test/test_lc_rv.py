import ellc
import numpy as np
import matplotlib.pyplot as plt


def main():
    period = 24.6
    t_obs = np.linspace(0, period, 300, endpoint=True)
    pars = [0.037, 0.061, 0.498, 91.642, 0, 19.592, 24.592, 47.9, 1.04, -0.15, 0.4, 0.38, 0.5, 0.32, 0.32]
    pars_rv = [0.037, 0.061, 0.498, 91.642, 0, 24.592, 47.9, 1.04, -0.15, 0.4, 0.38, 0.5, 0.32, 0.32]
    lc = ellc.lc(t_obs, *pars)
    rv1, rv2 = ellc.rv(t_obs, *pars_rv)
    plt.plot(t_obs, lc)
    plt.show()
    plt.plot(t_obs, rv1)
    plt.plot(t_obs, rv2)
    plt.show()


if __name__ == '__main__':
    main()