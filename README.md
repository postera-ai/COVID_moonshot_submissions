# Submissions to COVID Moonshot after second round

This repo contains all of the data submitted to the [COVID Moonshot](https://covid.postera.ai/covid) project as of Thursday, April 2nd, 11:59 PST, which marked the end of the second round of submissions.

The data has been split into 5 different files, in order to ease the triaging of the large number of compounds for synthesis. 

[covid_submissions_all_info.csv](https://github.com/mc-robinson/COVID_moonshot_submissions/blob/master/covid_submissions_all_info.csv)

This is the "master" file containing all of the submitted molecules. It contains quite a few fields described below:
- **SMILES**: the SMILES string as canonicalized by RDKit
- **CID**: the compound ID, with the first 3 strings corresponding to the submission ID at https://covid.postera.ai/covid/submissions. 
       e.g. "MAT-POS-ab1-3" is the 3rd molecule of submission "MAT-POS-ab1". The submission ID's are constructed as "FIRST 3 LETTERS OF            NAME - FIRST 3 LETTERS OF INSTITUTION - RANDOM 3 character string"
- **creator**: the creator of the submission
- **rationale**: the reasoning supplied by the creator for the molecule designs in the associated submission
- **fragments**: the fragments cited as inspiration for the submissions
- **real_space**: If the molecule is present in Enamine REAL space (note this is the 1.3 Billion molecule version, not the 13 Billion                         molecule version searchable through BioSolveIT
- **SCR**: If the molecule is in Enamine's Screening library
- **BB**: If the molecule is in Enamine's building blocks
- **link**: the link to the associated submission page
- **covalent_frag**: if one of the **fragments** cited as inpiration is a covalent fragment
- **covalent_warhead**: if *at least one* covalent warhead moiety is present in the designed molecule.
- **acrylamide**: If the acrylamide moiety is in the molecule
- **acrylamide_adduct**: If the acrylamide adduct moiety is in the molecule
- **chloroacetamide**: If the chloroacetamide moiety is in the molecule
- **chloroacetamide_adduct**: If the chloroacetamide adduct moiety is in the molecule 
- **vinylsulfonamide**: If the vinylsulfonamide moiety is in the molecule
- **vinylsulfonamide_adduct**: If the vinylsulfonamide adduct moiety is in the molecule
- **in_molport_or_mcule**: If the molecule orderable through Molport or Mcule, and not through Enamine, and thus is contained in the molport_or_mcule_df referenced below.
- **in_emolecules**: If the molecule is orderable through emolecules, but not through Enamine, Molport, or Mcule. It will thus be in the emolecules_df referenced below.

[enamine_purchaseable_df.csv](https://github.com/mc-robinson/COVID_moonshot_submissions/blob/master/enamine_purchaseable_df.csv)

This file contains all of the molecules that can be ordered or synthesized at Enamine. The relevant columns indicate the ID of the molecule if it is in REAL space, the Screening library, or Enamine's building blocks.

[molport_and_mcule_df.csv](https://github.com/mc-robinson/COVID_moonshot_submissions/blob/master/molport_and_mcule_df.csv)

This file contains all of the molecules *NOT available through Enamine* but purchasable through either Molport or Mcule.

[emolecules_df.csv](https://github.com/mc-robinson/COVID_moonshot_submissions/blob/master/emolecules_df.csv)

This file contains all of the molecule *NOT available through Enamine NOR Molport or Mcule*, but findable through emolecules.

[patent_df.csv](https://github.com/mc-robinson/COVID_moonshot_submissions/blob/master/patent_df.csv)

This file contains all of the molecules which cannot be directly purchased, but have been a part of reactions previously described in patents.








