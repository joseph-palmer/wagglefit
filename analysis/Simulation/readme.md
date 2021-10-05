## Honeybee foraging simulation
### Joseph Palmer <Joe.Palmer.2019@live.rhul.ac.uk>
#### Royal Holloway University of London

---

The simulation is a sequential agent based model (ABM).

An environment is drawn with flowers scattered randomly around a colony with a given intensity. Each flower has a resource level, location and distance to the colony.

Honeybees are composed of recruits and scouts. Scouts forage by moving from the colony in a linear path and reporting back the first flower in the path, a process equivilent to a 1D Poisson point process. Recruits select a flower resource provided by scouts and employed recruits based on a linear probability with distance of the resource to the colony, such that closer resources have a higher probability of being selected than those futher away. One they have a resource the recruits will forage at the source until the flower is depleted (resourse level goes to 0) or until the flower is randomly removed from the environment. After foraging the recruits add their flower location to the pool of dances unemployed recruits sample from.

Scouts and recruits forage in each iteration and the distances danced are recorded for analysis. Each iteration $n$ flowers are randomly removed from the environment and $k$ flowers are added.
