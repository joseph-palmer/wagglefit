#!/usr/bin/env python3
"""Classes of agents (workers, scouts, flowers)"""
__appname__  = "Agents.py"
__author__   = "Joseph Palmer <joseph.palmer18@imperial.ac.uk>"
__version__  = "0.0.1"
__date__     = "12-2019"

import numpy as np
import matplotlib.pyplot as plt
import Methods
import config
import copy
from shapely.geometry import Point, Polygon
import shapely

class Flower:
    # Attributes #
    location = (None, None)
    location_point = None
    distance = None
    resource_level = None
    exists = 0

    # Methods #
    def __init__(self, location, quality_array):
        self.location = location
        self.location_point = Point(*location)
        self.distance = Methods.getdist(*location)
        self.resource_level = 2000 #np.random.randint(10, 5000)
        self.quality_array = quality_array #quality_array[np.random.randint(0, len(quality_array))]
        self.exists = 1
        return None

    def remove_flower(self):
        self.location = (None, None)
        self.distance = None
        self.resource_level = None
        self.exists = 0
        return None

    def visit(self):
        if self.resource_level > 0:
            self.resource_level -= 1
        else:
            self.remove_flower()
        return None

class Honeybee:
    flower = None
    distance = None
    mode = None # -1 is scout, 1 is employed recruit, 0 is unemployed recruit

class Scout(Honeybee):
    def __init__(self):
        self.mode = -1
        return None

    def _drawpath(self, s1, s2, flowers, infa=[], i=0, limmit=10, plot=True):
        angle = np.random.uniform(0, 360)
        distance = np.random.uniform(0, 1)
        base_fa = [
            tuple(s1),
            (s1[0], distance+s1[1]),
            (s2[0], distance+s2[1]),
            tuple(s2)
        ]
        fa = Polygon(base_fa)
        fa = shapely.affinity.rotate(fa, angle, origin=s1)
        if plot:
            x,y = fa.exterior.xy
            plt.plot(x,y, color="red")
        a, b = list(fa.exterior.coords)[1:3]
        infa_idx = [fa.contains(j.location_point) for j in flowers]
        [infa.append(ii) for ii in np.array([j.distance for j in flowers])[infa_idx]]
        if i < limmit:
            self._drawpath(a, b, flowers=flowers, infa=infa, i=i+1, plot=plot)
        return infa

    def search(self, plot=False):
        """
        Scout leaves hive in given direction covering a given area. Retrive
        all flowers in thar area and select the one closest to the hive.

        Arguments:
            HoneyBee {Honeybee.object} -- The honeybee
            plot {bool} -- Adds the search space to the plot (default False)
        """
        # extract flowers if they are in the environment
        #idex = list(config.FLOWERDICT.keys())
        flowers = list(config.FLOWERDICT.values())
        idex = [i for i in range(len(flowers)) if flowers[i].exists]
        flowers = np.array(flowers)[idex]

        # random walk scout searching
        infa = self._drawpath(s1=[0.01, 0.01], s2=[-0.01, -0.01], flowers=flowers, infa=[], i=0)
        if len(infa) > 0:
            # get the minimum value in foraging area and visit flower
            vfunc = np.vectorize(lambda x: x.distance)
            match = np.where(vfunc(flowers) == min(infa))[0][0]
            self.distance = min(infa)
            self.flower = flowers[match]
            self.flower.visit()
            self.update()
        else:
            self.flower = None
            self.distance = None

        if plot:
            if self.flower:
               plt.plot(*self.flower.location, "o", color="g")

        return None

    def update(self):
        if self.flower.exists:
            pass
        else:
            self.flower = None
            self.distance = None
        return None

class Recruit(Honeybee):
    def __init__(self):
        self.mode = 0
        return None

    def forage(self):
        if self.flower.exists:
            self.flower.visit()
        else:
            self.update(None, None)
        return None

    def update(self, flower, distance):
        if flower and distance:
            self.mode = 1
            self.flower = flower
            self.distance = distance
        else:
            self.mode = 0
            self.flower = None
            self.distance = None
        return None

    def watch_dance(self):
        keys = list(config.BROADCASTED_DANCES.keys())
        vals = list(config.BROADCASTED_DANCES.values())
        if len(vals) > 0:
            idx = np.argsort(vals)
            probs = np.linspace(1, 0, len(vals))
            probs = probs/sum(probs)
            mflower = np.random.choice(np.array(keys)[idx], p=probs)
            mdist = config.BROADCASTED_DANCES[mflower]
            self.update(mflower, mdist)
        return None

    def sample_dance_random(self):
        # select randomly from all dances
        flowers = np.array(list(config.BROADCASTED_DANCES.keys()))
        if len(flowers) > 0:
            mflower = np.random.choice(flowers)
            mdist = config.BROADCASTED_DANCES[mflower][0]
            self.update(mflower, mdist)
            self.forage()

    def sample_dance_bias_dq(self):
        # sample according to quality
        flowers = np.array(list(config.BROADCASTED_DANCES.keys()))
        if len(flowers) > 0:
            distances = np.array(list(config.BROADCASTED_DANCES.values()))[:,0]
            quality = np.array(list(config.BROADCASTED_DANCES.values()))[:,-1]
            qd = Methods.qual2dance(quality, distances, 10)
            p = qd/sum(qd)
            idx = np.arange(0, len(flowers))
            if not np.isnan(np.sum(p)):
                mflower = flowers[int(np.random.choice(idx.astype(float), p=p.astype(float)))]
                mdist = config.BROADCASTED_DANCES[mflower][0]
                self.update(mflower, mdist)
                self.forage()

    def sample_dance_best(self):
        # sample according to quality
        flowers = np.array(list(config.BROADCASTED_DANCES.keys()))
        if len(flowers) > 0:
            distances = np.array(list(config.BROADCASTED_DANCES.values()))[:,0]
            quality = np.array(list(config.BROADCASTED_DANCES.values()))[:,-1]
            qd = Methods.qual2dance(quality, distances, 10)
            bestidx = np.argmax(qd)
            mflower = flowers[bestidx]
            mdist = config.BROADCASTED_DANCES[mflower][0]
            self.update(mflower, mdist)
            self.forage()


    # def search(self, plot=False):
    #     """
    #     Scout leaves hive in given direction covering a given area. Retrive
    #     all flowers in thar area and select the one closest to the hive.

    #     Arguments:
    #         HoneyBee {Honeybee.object} -- The honeybee
    #         plot {bool} -- Adds the search space to the plot (default False)
    #     """
    #     # extract flowers if they are in the environment
    #     #idex = list(config.FLOWERDICT.keys())
    #     flowers = list(config.FLOWERDICT.values())
    #     idex = [i for i in range(len(flowers)) if flowers[i].exists]
    #     flowers = np.array(flowers)[idex]

    #     # place a rectangle from the hive outwards of a given width in a random
    #     # direction, get all flowers in the rectangle (foraging area)
    #     angle = np.random.uniform(0, 360)
    #     base_fa = [
    #         (0.01, 0.01),
    #         (0.01, 1.01),
    #         (-0.01, 1.01),
    #         (-0.01, -0.01)
    #     ]

    #     fa = Polygon(base_fa)
    #     fa = shapely.affinity.rotate(fa, angle, origin=(0,0))
    #     infa_idx = [fa.contains(i.location_point) for i in flowers]
    #     infa = np.array([i.distance for i in flowers])[infa_idx]
    #     if len(infa) > 0:
    #         # get the minimum value in foraging area and visit flower
    #         vfunc = np.vectorize(lambda x: x.distance)
    #         match = np.where(vfunc(flowers) == min(infa))[0][0]
    #         self.distance = min(infa)
    #         self.flower = flowers[match]
    #         self.flower.visit()
    #         self.update()
    #     else:
    #         self.flower = None
    #         self.distance = None

    #     # show point in foraging area
    #     if plot:
    #         x,y = fa.exterior.xy
    #         endp1, endp2 = list(fa.exterior.coords)[1:3]
    #         self._drawpath(endp1, endp2)
    #         plt.plot(x,y)
    #         if self.flower:
    #             plt.plot(*self.flower.location, "o", color="g")

    #     return None
