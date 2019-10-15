# genfusgen

Construction of Generalized Fusenes

## Usage

```
genfusgen [-p] [-d] [-k] [-f] [-o OUTFILE] [-m M] [-i I] SPECS

 -p,--planarcode write planar code to stdout or outfile
 -d,--duals      generate inner duals
 -r,--regular    filter subgraphs of regular lattice
 -k,--kekule     filter kekule structures
 -f,--fix        fix the outer face in clockwise direction of edge (1,2)
 -o,--output     write to OUTFILE instead of stdout
 -m,--modulo
 -i,--index      only use inner dual with index I (modulo M)
 SPECS           sequence of pairs "n:m", meaning there are m n-gons
```
