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

### Tutorial

This will walk through the design of a 6 helix bundle triangle with continuous, rounded corners.

First, of all, we need to know our exact dimensions before creating a template.

In this design, we'll be folding a 6 helix bundle along its 3-plane axis. We can mark the planes as such:

<img src="/EXAMPLES/walkthrough/6hbplanes.png" alt="6 helix bundle planes" width="50%" />

We then have to decide how much we want our corners to curve.

![](/EXAMPLES/walkthrough/curvedtrianglediagram.png)

In this case, we are only interested in the inner dimensions of the triangle, but the other layers have been drawn with dotted lines. The longer dashed line at the triangle's corners are how much we will be shortening an edge to insert the curved corner. We started with a 126 bp edge.

For this structure, with a bit of trig, we determined that, along our bend, our inner layer should have 32 bp (approximately 3 full turns).

Using the beam model equation above, for the positions of our three planes, the insertion/deletion gradient is [-14, 0, +14]. Thus, with 32 bp as our reference, the subsequently outer layer bends are 46 bp and 60 bp long.

We use the middle plane length as the nominal length for creating the caDNAno template. The total length of our triangle (as a cycle) is $46*3+84*3 = 390$.

We then go on and draw that in caDNAno.

![](/EXAMPLES/walkthrough/step1.png)

Add our breakpoints...

![](/EXAMPLES/walkthrough/step2.png)

Note that the endpoints of the structure are connected.

<img src="/EXAMPLES/walkthrough/edge1.png" alt="Edge 1" width="50%" /><img src="/EXAMPLES/walkthrough/edge2.png" alt="Edge 2" width="50%" />

We now run the script, three times, once to make each corner.

`python autobender.py -o temp_01.json -lt hc -a 120 -l 46 -s 91 -x 0 template.json `

Produces...

![](/EXAMPLES/walkthrough/temp_01.png)

For the second bend...

`python autobender.py -o temp_02.json -lt hc -a 120 -l 46 -s 221 -x 0 temp_01.json `

![](/EXAMPLES/walkthrough/temp_02.png)

And for the third bend...

`python autobender.py -o temp_03.json -lt hc -a 120 -l 46 -s 351 -x 0 temp_02.json `

![](/EXAMPLES/walkthrough/temp_03.png)

This final caDNAno design can then be converted by [tacoxDNA](http://tacoxdna.sissa.it/) and simulated in [oxDNA](https://dna.physics.ox.ac.uk/index.php/Main_Page) or viewed in [oxView](https://sulcgroup.github.io/oxdna-viewer/).

![](/EXAMPLES/walkthrough/canvas.png)

### Links
CanDo: https://cando-dna-origami.org/

### TODO: 
- ~Add an input parameter for rotating the desired normal plane axis.~ Added!
- Flag some known structural issues (overly dense insertion/deletion gradient)

### References
[^1]: Dietz, H., Douglas, S. M., & Shih, W. M. (2009). Folding DNA into twisted and curved nanoscale shapes. Science, 325(5941), 725-730. [[Science](https://www.science.org/doi/10.1126/science.1174251)]
