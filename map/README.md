# RBC-GEM map files

This directory contains the RBC-GEM map files.

## Recommended Use
As stated in the Escher [documentation](https://escher.readthedocs.io/en/latest/tips-and-tricks.html#escher-performance-with-large-maps):

>_Escher performance with large maps
Escher works best with maps that have less than about 200 reactions. If you are working with more reactions, we recommend splitting your map into multiple small maps. A trick for splitting up a map is to first select the reactions you want to keep, then choose “Edit > Invert Selection” and then “Edit > Delete”._

As such, recommended use of the RBC map is therefore to generate smaller maps of the subsystems of interest, then add the corresponding transporters and exchanges.

The following can help improve performance should it be desiredable to use the larger map

>* _Turn off tooltips in the settings menu, especially tooltips over Objects._
>* _Try the “Use 3D Transform” option in the settings menu. Depending on your browser, this can increase responsiveness when moving around the map._
>* _Turn off labels, gene reaction rules, and/or secondary metabolites in the settings menu._

See the Escher documentation for more information.

## File Descriptions

* _RBC-GEM.full.map.json_: The map JSON file created by Escher.
* _RBC-GEM.full.map.xml_: The map as an SBML layout file of the map, created by the [EscherConverter software](https://escher.readthedocs.io/en/latest/escherconverter.html#escherconverter).
* _RBC-GEM.full.map.sbgn_: The map as an SBGN file of the map, created by the [EscherConverter software](https://escher.readthedocs.io/en/latest/escherconverter.html#escherconverter).
* _RBC-GEM.full.map.pdf_: A PDF version of the full map with reactions color coded based on general metabolic category. Color codes correspond to Figure 2 in the manuscript. List Numbers also correspond to categories used.
    1. Amino Acids (yellow)
    2. Carbohydrates (dark green)
    3. Lipids (blue)
    4. Cofactors and Vitamins (light green)
    5. Nucelotides (aquamarine)
    6. Reactive species (red)
    7. Other (black)
    8. Transport reactions (purple)
