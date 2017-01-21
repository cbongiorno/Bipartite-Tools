# Bipartite-Tools
A tool for testing community detection algorithms on one mode bipartite projected network.

To install on Linux run:

```
sudo ./INSTALL
```

Then It is possible to use as a python library:

```
import bipnet
```

The module validate contains the functions to estimate the statistically validated networks, the module metrics contains other useful functions. The code includes a wrapper to the louvain community detection method written by E. Lefebvre and released under GNU Licence. The louvain code is available separately [here](https://sourceforge.net/projects/louvain/)
