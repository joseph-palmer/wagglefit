#!/usr/bin/env python3
"""Run simulation of honeybee foraging"""
__appname__  = "run.py"
__author__   = "Joseph Palmer <joseph.palmer18@imperial.ac.uk>"
__version__  = "0.0.2"
__date__     = "02-2020"

import numpy as np
import pandas as pd
import scipy.stats
import matplotlib.pyplot as plt
import Agents
import Methods
import config

def create_flowers(lwr, upr, res, quality_array):
    # Simulation window parameters
    r = upr-lwr
    xx0 = 0
    yy0 = 0

    areaTotal = np.pi*r**2

    # Point process parameters
    lambda0 = res

    # Simulate Poisson point process
    numbPoints = scipy.stats.poisson(lambda0*areaTotal).rvs()
    theta = 2*np.pi*scipy.stats.uniform.rvs(0,1,((numbPoints,1)))
    rho = r*np.sqrt(scipy.stats.uniform.rvs(0,1,((numbPoints,1))))

    # Convert from polar to Cartesian coordinates
    xx = rho * np.cos(theta)
    yy = rho * np.sin(theta)

    # Shift centre of disk to (xx0,yy0)
    xx=xx+xx0
    yy=yy+yy0

    # Set locations as flowers
    locations = np.column_stack((xx, yy))
    flowers = np.apply_along_axis(Agents.Flower, 1, locations, quality_array)
    fdict = {
        i : flowers[i] for i in range(res)
    }
    return fdict

def initialize(lwr, upr, res, nbees, nscouts, quality_array):
    # create flowers in environment and store globaly.
    config.FLOWERDICT = create_flowers(lwr, upr, res, quality_array)

    # create honeybee population
    nrecruits = nbees-nscouts
    scouts = np.array([Agents.Scout() for i in range(nscouts)])
    recruits = np.array([Agents.Recruit() for i in range(nrecruits)])
    return (scouts, recruits)

def flower_turnover(n, lwr, upr, k, quality_array):
    # select n flowers and remove
    assert (n <= len(config.FLOWERDICT.keys())), "N > flowers in environment!"
    rmv_idx = np.random.choice(list(config.FLOWERDICT.keys()), n, replace=False)
    [config.FLOWERDICT[i].remove_flower() for i in rmv_idx]

    # add k new flowers
    newflower = create_flowers(lwr, upr, k, quality_array)
    config.FLOWERDICT.update(newflower)
    return None

def simulate(lwr, upr, res, nbees, nscouts, nsim, tn, tk, quality_array, plot=False, sampletype=0):
    # set up environment
    scouts, recruits = initialize(lwr, upr, res, nbees, nscouts, quality_array)
    recruit_distances = []
    scout_distances = []

    # set up for plot info
    fx = []
    fy = []
    dist = []
    sx = []
    sy = []
    rx = []
    ry = []


    for j in range(nsim):
        flower_turnover(tn, lwr, upr, tk, quality_array)

        # plot environment
        flowers = list(config.FLOWERDICT.values())
        idex = [i for i in range(len(flowers)) if flowers[i].exists]
        flowers = np.array(flowers)[idex]

        # plot enviroment
        if plot:
            x = [i.location[0] for i in flowers]
            y = [i.location[-1] for i in flowers]
            plt.plot(x, y, "o", color="b")

        if len(flowers) < 1:
            exit("No flowers left to forage at.")

        # scouts search, employed recruits forage
        [i.search(plot) for i in scouts]
        spoints = [j.flower.location for j in scouts if j.flower and j.flower.resource_level]
        rpoints = [j.flower.location for j in recruits if j.mode > 0]

        # record flower info
        [(fx.append(i.location[0]), fy.append(i.location[-1])) for i in flowers]
        [(dist.append(i.distance)) for i in flowers]
        [(sx.append(j[0]), sy.append(j[-1])) for j in spoints]
        [(rx.append(j[0]), ry.append(j[-1])) for j in rpoints]

        # recruitment
        s_dist = [i.distance for i in scouts if i.flower and i.flower.resource_level]

        r_dist = [i.distance for i in recruits if i.flower and i.flower.resource_level]
        dists = s_dist + r_dist
        s_qual = [i.flower.quality_array for i in scouts if i.flower and i.flower.resource_level]
        r_qual = [i.flower.quality_array for i in recruits if i.flower and i.flower.resource_level]
        quals = s_qual + r_qual
        s_dance = [i.flower for i in scouts if i.flower  and i.flower.resource_level]
        r_dance = [i.flower for i in recruits if i.flower and i.flower.resource_level]
        dances = s_dance + r_dance

        config.BROADCASTED_DANCES = {
            dances[i] : [dists[i], quals[i]] for i in range(len(dances))
        }

        # unemployed watch the dances
        if sampletype == 0:
            [i.sample_dance_random() for i in recruits]
        elif sampletype == 1:
            [i.sample_dance_bias_dq() for i in recruits]
        else:
            [i.sample_dance_best() for i in recruits]

        # update
        recruit_distances += r_dist
        scout_distances += s_dist

    return (scout_distances, recruit_distances, [fx, fy, sx, sy, rx, ry, dist])


def main():
    # paramaters for model
    nruns = 5
    lwr = -2.5
    upr = 2.5
    res = 5000
    nbees = 100
    nscouts = 20
    nsim = 100
    tn = 10
    tk = 20
    quality_array = np.random.uniform(0, 10, 10)


    # run simulation
    print(f"Running simulation for {nsim} iterations")
    total_sdist = np.array([])
    total_rdist = np.array([])
    for i in range(nruns):
        print(f"\t - {i+1}/{nruns}")
        sdist, rdist, additional_params = simulate(
            lwr,
            upr,
            res,
            nbees,
            nscouts,
            nsim,
            tn,
            tk,
            quality_array,
            sampletype=1,
            plot=False
        )
        total_rdist = np.concatenate(
            [total_rdist, np.random.choice(rdist, 60, replace = False)]
        )
        total_sdist = np.concatenate(
            [total_sdist, np.random.choice(sdist, 60, replace = False)]
        )

    # sdist_path = "../Results/scout_distribution_v2.csv"
    # rdist_path = "../Results/recruit_distribution_v2.csv"
    # dists = [total_sdist, total_rdist]
    # paths = [sdist_path, rdist_path]
    # for i in range(2):
    #     np.savetxt(paths[i], dists[i], delimiter = ",")
    return 0

if __name__ == "__main__":
    main()
