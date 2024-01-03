# RBC-GEM data

This directory contains the RBC-GEM data files.

### databases
Contains useful files for mapping model content to different databases. Files stored here are intended to be used repeatedly and updated with the corresponding databases.

### deprecatedIdentifiers
Contains the `tsv` files with reaction, metabolite, and gene annotations, identical to that in `/model/` for those identifiers that have been used in previous iterations RBC-GEM and subsequently removed. The identifiers are kept here so as to avoid accidentally reusing them in the future.

### external
Contains `.zip` files for various publications and unmodified supplemental files assoicated with datasets used in working with the RBC-GEM reconstruction. Even publications and their associated data may be unavailable over time (e.g., change to URL for a database). Archiving the relevant publications and publically available helps alleviate these difficulties. Intended only for academic purposes.

### memote
Contains data relevant to MEMOTE experiments

### raw
Contains the raw RBC-GEM data files. It is intended to be used to store raw/initial data files that will be utilized, transformed, and/or processed in some capacity using subsequent code. It should remain empty on the `main` branch. 

### interim
Contains the interim RBC-GEM data files. It is intended to be used to store data that is being processed/finalized. It should remain empty on the `main` branch. 

### processed
Contains the processed RBC-GEM data files. It is intended to be used to store processed/finalized data before organization into the appropriate repository locations. It should remain empty on the `main` branch. 