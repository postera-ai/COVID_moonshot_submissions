# COVID Moonshot Data

[![License: CC0-1.0](https://licensebuttons.net/l/zero/1.0/80x15.png)](http://creativecommons.org/publicdomain/zero/1.0/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repo contains all of the designs submitted to the [COVID Moonshot](https://covid.postera.ai/covid) project, as well as the initial experimental data created as part of the project.

The data has been split into many different folders and files, in order to ease the organization of such a large number of designs made through several CROs.

There exists one "master file" containing most of the non-experimental information: [covid_submissions_all_info.csv](https://github.com/mc-robinson/COVID_moonshot_submissions/blob/master/covid_submissions_all_info.csv) This file is described in much greater detail at the end of this README.

There also exist several folders:
- [data_for_CDD/](https://github.com/mc-robinson/COVID_moonshot_submissions/tree/master/data_for_CDD): prepares the data in a format suitable for upload to the CDD (Collaborative Drug Discovery) vault
- [orders/](https://github.com/mc-robinson/COVID_moonshot_submissions): The compounds that have been ordered up to this point.
- [shipments/](https://github.com/mc-robinson/COVID_moonshot_submissions/tree/master/shipments_data) the compounds that have been recevied up to this point.

Additionally, the remote [Moonshot_DR_curves](https://drive.google.com/drive/folders/1qhhDSImiu2f-5IiI0sI5OfSTxeM8JQ2A?usp=share_link) folder contains the biochemical dose-response curves for the compounds on this project.

Some more detail on [covid_submissions_all_info.csv](https://github.com/mc-robinson/COVID_moonshot_submissions/blob/master/covid_submissions_all_info.csv), and its many fields.

This is the main file containing all of the submitted molecules. It contains quite a few fields described below:
- **SMILES**: the SMILES string as canonicalized by RDKit and standardized by https://github.com/chembl/ChEMBL_Structure_Pipeline. Note that this SMILES is the registration SMILES, but may not correspond to the actual SMILES of the compound. See the `suspected_SMILES` field of `suspected_SMILES.csv` for further info.
- **CID**: the compound ID, with the first 3 strings corresponding to the submission ID at https://covid.postera.ai/covid/submissions. 
       e.g. "MAT-POS-916a2c5a-1" is the 3rd molecule of submission "MAT-POS-916a2c5a". 
       The submission ID's are constructed as "<FIRST 3 LETTERS OF NAME>-<FIRST 3 LETTERS OF INSTITUTION>-<RANDOM 8 character string>"
- **CID (canonical)**: The canonical compound ID corresponding to the ID given when the       molecule was submitted for the first time. 
- **CID (old format)**: The CID in the old format in order to maintain backwards compatability. This format bas been deprectated.
- **PostEra Link**: the link to the associated molecule detail page on postera.ai
- **InChIKey**: The InChIKey of the molecule.
- **creator**: the creator of the submission. 
- **rationale**: the reasoning supplied by the creator for the molecule designs in the 
associated submission
- **Submission Notes**: Any notes supplied by the creator with the submission.
- **fragments**: the XChem fragments cited as inspiration for the submissions. This field has now been deprecated for the "Inspired By" field found on the website.
- **Structure ID**: The Fragalysis ID of the crystal structure, if it exists. 
- **Fragalysis Link**: The link to the structure on Fragalysis, if it exists
- **MW**: Molecular weight of the molecule
- **cLogP**: Calculated LogP according to RDKit
- **HBD**: number of hydrogen bond donors according to RDKit
- **HBA**: number of hydrogent bond acceptors according to RDKit
- **TPSA**: Topological polar surface area according to RDKit
- **series**: The "chemical series of interest" we consider the molecule to be a part of.
- **CDD Name**: The CDD name displayed in the CDD vault.
- **CDD_mol_ID**: The molecule ID in CDD, useful mostly for the API.
- **r_inhibition_at_20_uM**: Single shot percent inhibition of molecule at 20 uM in RapidFire assay
- **r_inhibition_at_50_uM**: Single shot percent inhibition of molecule at 50 uM in RapidFire assay
- **r_avg_IC50**: Average IC50 over the multiple dose-response runs of the molecule in the Rapidfire assay
- **r_curve_IC50**: List of IC50 values from each dose-response run of molecule in the RapidFire assay
- **r_max_inhibition_reading**: List of the the maximum inhibition value from each dose-response run of the RapidFire assay.
- **r_min_inhibition_reading**: List of the the minimum inhibition value from each dose-response run of the RapidFire assay.
- **r_hill_slope**: List of the the Hill slopes from the dose-response curves for each run of the RapidFire assay.
- **r_R2**: List of the the R-squared values from the dose-response curves for each run of the RapidFire assay.
- **r_concentration_uM**: List of lists. Each inner list consists of the concentrations tested in the RapidFire assay for dose-response.
- **r_inhibition_list**: List of lists. Each inner list consists of the inhibition values at the concentrations in r_concentration_uM
- **f_inhibition_at_20_uM**: Single shot percent inhibition of molecule at 20 uM in Fluorescence assay
- **f_inhibition_at_50_uM**: Single shot percent inhibition of molecule at 50 uM in Fluorescence assay
- **f_avg_IC50**: Average IC50 over the multiple dose-response runs of the molecule in the Fluorescence assay
- **f_avg_pIC50**: Average pIC50 over the multiple dose-response runs of the molecule in the Fluorescence assay
- **f_curve_IC50**: List of IC50 values from each dose-response run of molecule in the Flourescence assay
- **f_max_inhibition_reading**: List of the the maximum inhibition value from each dose-response run of the Fluorescence assay.
- **f_min_inhibition_reading**: List of the the minimum inhibition value from each dose-response run of the Fluorescence assay.
- **f_hill_slope**: List of the the Hill slopes from the dose-response curves for each run of the Fluorescence assay.
- **f_R2**: List of the the R-squared values from the dose-response curves for each run of the Fluorescence assay.
- **f_concentration_uM**: List of lists. Each inner list consists of the concentrations tested in the Fluorescence assay for dose-response.
- **f_inhibition_list**: List of lists. Each inner list consists of the inhibition values at the concentrations in f_concentration_uM
- **relative_solubility_at_20_uM**: Nephelometry-based solubility assay to define threshold compound solubility at 20 uM
- **relative_solubility_at_100_uM**: Nephelometry-based solubility assay to define threshold compound solubility at 100 uM
- **trypsin_IC50**: IC50 in the Trypsin counter-assay.
- **NMR_std_ratio**: Saturation Transfer Difference (STD) NMR experiments of Mpro with ligands. Higher ratio indicates binding.
- **ORDERED**: If the compound has been ordered for synthesis.
- **ORDER_DATE**: The date on which the molecule was ordered for synthesis / shipment.
- **MAKER**: Where the molecule was ordered from on the Order Date.
- **MADE**: If the ordered compound has been synthesized and delivered.
- **SHIPMENT_DATE**: The date on which the molecule was shipped to us.
- **ASSAYED**: If the molecule has been tested for inhibition.

## Details on SMILES and molecule identities##

Please see `suspected_SMILES.csv`: Note that many of the enantiopure compounds on this project were obtained by chiral separation, and thus the compounds are often obtained as a single enantiomer with unknown absolute stereochemistry. The `suspected_SMILES` column represents the understanding of the actually identity of the compound given the current information, which explains the frequent use of enhanced/relative stereochemistry representations. Keep in mind that the`suspected_SMILES` is often not the same as the `registration_SMILES`, since compounds must be registered as a real (non-relatively defined) compound on this project.

## License

All code is [MIT licensed](LICENSE) and all data is [CC0 1.0](https://creativecommons.org/publicdomain/zero/1.0/)
