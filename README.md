# Bipartite-Tools
A tool for testing community detection algorithms on bipartite projected networks.

To install on Linux run:

```
sudo ./INSTALL
```

It is possible to use it as a python library:

```
import bipnet
```

The module validate contains the functions for statistically validated networks. The module metrics contains other useful functions (for ex. the AWI). The code includes a wrapper to the louvain community detection method written by E. Lefebvre and released under GNU Licence. The louvain code is available separately [here](https://sourceforge.net/projects/louvain/)
