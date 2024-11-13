# dnahbbender

Includes some simple scripts to add insertions and deletions automatically for creating bent DNA origami helix bundles. Only works with caddnano2's JSON file format. Does not work with the data structure of cadnano2.5.

This is a geometric model only. It does not implement the energy minimization model from [^1].

Insertion/deletion gradient is calculated using a beam model that simply implements:

$\Delta bp = \theta \cdot d \cdot \delta_i$

  - $\Delta bp$: # of insertions or deletions
  - $\theta$: Desired bend angle (degrees)
  - $d$: Axial rise (DEFAULT = 0.332 nm)
  - $\delta_i$: Distance from normal plane

Supports both honeycomb (`hc`) and square (`sq`) lattice designs.

Does not need to be compiled. Download and run the script using Python:

`python autobender.py -f <input filename> -o <output filename> -lt <lattice type> -a <bend angle> -l <bend region length> -s <starting index>`

e.g., `python autobender.py -f 6helixbundle.json -o output.json -lt hc -a 120 -l 90 -s 211` will access the 6helixbundle.json file designed using the honeycomb lattice, start adding the insertions/deletions at cursor index 211 (the yellow bar in cadnano), in a region of 90 bps (211 to 300), and enough to create a 120° bend angle, then save the new JSON data in output.json.

Insertions and deletions are distributed uniformly. A correction script is also ran to shift edits off crossover positions (which will cause problems and the CanDo simulation will fail).

:point_right:**Note**, the script cannot discern between whether the input file is actually of a honeycomb or square lattice. The script will still appear to work, but the number of insertions/deletions added will not produce the desired bend angle.

The repository also includes an assortment of previously created or processed templates and examples. Structures primarily verified with CanDo.

:point_right:**Outputs of this model have not been verified experimentally.**

### Links
CanDo: https://cando-dna-origami.org/

### TODO: 
- Add an input parameter for rotating the desired normal plane axis.
- Flag some known structural issues (overly dense insertion/deletion gradient)

### References
[^1]: Dietz, H., Douglas, S. M., & Shih, W. M. (2009). Folding DNA into twisted and curved nanoscale shapes. Science, 325(5941), 725-730. [[Science](https://www.science.org/doi/10.1126/science.1174251)]
