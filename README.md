# COVID Moonshot Data

This repo contains all of the data submitted to the [COVID Moonshot](https://covid.postera.ai/covid) project.
The data has been split into many different folders and files, in order to ease the triaging of the large number of compounds for synthesis.

Each folder will contain its own README file.  PLEASE READ THESE if you wish to contribute data to the project. NO ONE should be spending time standardizing multiple fields (SMILES/smiles/Smiles/smile_string). Please follow *the exact speifications* for all future files, or else a program will reject them. All future files must also include a SMILES field.

There exists one "master file" containing most of the non-experimental information: [covid_submissions_all_info.csv](https://github.com/mc-robinson/COVID_moonshot_submissions/blob/master/covid_submissions_all_info.csv) This file is described in much greater detail at the end of this README.

There also exist several folders:
- [availability data/](https://github.com/mc-robinson/COVID_moonshot_submissions/tree/master/availability_data): contains files identifying the molecules that are avaiable from Enamine, Molport and/or Mcule, and eMolecules.
- [data_for_CDD/](https://github.com/mc-robinson/COVID_moonshot_submissions/tree/master/data_for_CDD): prepares the data in a format suitable for upload to the CDD (Collaborative Drug Discovery) vault
- [experimental_data/](https://github.com/mc-robinson/COVID_moonshot_submissions/tree/master/experimental_data): data resulting from assays on the ordered compounds.
- [orders_data/](https://github.com/mc-robinson/COVID_moonshot_submissions): The compounds that have been ordered up to this point.
- [shipments_data/](https://github.com/mc-robinson/COVID_moonshot_submissions/tree/master/shipments_data) the compounds that have been recevied up to this point.
- [submissions_data/](https://github.com/mc-robinson/COVID_moonshot_submissions/tree/master/submissions_data) specialiazed analysis of all info that has been submitted to covid.postera.ai/covid. 

Some more detail on [covid_submissions_all_info.csv](https://github.com/mc-robinson/COVID_moonshot_submissions/blob/master/covid_submissions_all_info.csv), and its many fields.

This is the "master" file containing all of the submitted molecules. It contains quite a few fields described below:
- **SMILES**: the SMILES string as canonicalized by RDKit and standardized by https://github.com/chembl/ChEMBL_Structure_Pipeline.
- **CID**: the compound ID, with the first 3 strings corresponding to the submission ID at https://covid.postera.ai/covid/submissions. 
       e.g. "MAT-POS-ab1-3" is the 3rd molecule of submission "MAT-POS-ab1". The submission ID's are constructed as "FIRST 3 LETTERS OF            NAME - FIRST 3 LETTERS OF INSTITUTION - RANDOM 3 character string"
- **creator**: the creator of the submission
- **rationale**: the reasoning supplied by the creator for the molecule designs in the associated submission
- **fragments**: the fragments cited as inspiration for the submissions
- **link**: the link to the associated submission page
- **real_space**: If the molecule is present in Enamine REAL space (note this is the 1.3 Billion molecule version, not the 13 Billion                         molecule version searchable through BioSolveIT
- **SCR**: If the molecule is in Enamine's Screening library
- **BB**: If the molecule is in Enamine's building blocks
- **extended_real_space**: If the molecule is in the full 13 billion molecule real space
- **in_molport_or_mcule**: If the molecule orderable through Molport or Mcule, and not through Enamine, and thus is contained in the molport_or_mcule_df referenced below.
- **in_emolecules**: If the molecule is orderable through emolecules, but not through Enamine, Molport, or Mcule. It will thus be in the emolecules_df referenced below.
- **covalent_frag**: if one of the **fragments** cited as inpiration is a covalent fragment
- **covalent_warhead**: if *at least one* covalent warhead moiety is present in the designed molecule.
- **acrylamide**: If the acrylamide moiety is in the molecule
- **acrylamide_adduct**: If the acrylamide adduct moiety is in the molecule
- **chloroacetamide**: If the chloroacetamide moiety is in the molecule
- **chloroacetamide_adduct**: If the chloroacetamide adduct moiety is in the molecule 
- **vinylsulfonamide**: If the vinylsulfonamide moiety is in the molecule
- **vinylsulfonamide_adduct**: If the vinylsulfonamide adduct moiety is in the molecule
- **nitrile**: If the nitrle moiety is in the molecule
- **nitrile_adduct**: If the nitrle adduct moiety is in the molecule
- **MW**: Molecular weight of the molecule
- **cLogP**: Calculated LogP according to RDKit
- **HBD**: number of hydrogen bond donors according to RDKit
- **HBA**: number of hydrogent bond acceptors according to RDKit
- **TPSA**: Topological polar surface area according to RDKit
- **num_criterion_violations**: The number of violations of the rough criterion set out by the med-chem team:
```
{
    "MW": [0, 550],
    "cLogP": [-1, 5],
    "HBD": [0, 5],
    "HBA": [0, 10],
    "TPSA": [0, 200],
}
```
- **BMS**: If the molecule passes the BMS (Bristol-Myers Squibb) alerts described in DOI: 10.1021/ci050504m
- **Dundee**: If the molecule passes the Dundee alerts described in DOI: 10.1002/cmdc.200700139
- **Glaxo**: If the molecule passes the Glaxo-Wellcome alerts described in https://doi.org/10.1021/ci990423o
- **Inpharmatica**: If the molecule passes the Inpharmatica alerts described in https://www.surechembl.org/knowledgebase/169485
- **LINT**: If the molecule passes the Pfizer LINT alerts described in DOI: 10.2174/157340605774598081
- **MLSMR** If the molecule passes the MLSMR filters described in https://mlsmr.evotec.com/MLSMR_HomePage/pdf/MLSMR_Excluded_Functionality_Filters_200605121510.pdf
- **PAINS**: If the molecules passes the PAINS filters described in DOI: 10.1021/jm901137j
- **SureChEMBL**: If the molecule passes the SureChEMBL filters described in https://www.surechembl.org/knowledgebase/169485
- **ORDERED**: If the compound has been ordered for synthesis.
- **MADE**: If the ordered compound has been synthesized and delivered.

Special thanks goes to Pat Walters for compiling the SMARTS for the alerts here: https://github.com/PatWalters/rd_filters

The above master file is also subdivided into many files, which are described in the respective sub-folders.









