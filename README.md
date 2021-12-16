# COVID Moonshot Data

This repo contains all of the designs submitted to the [COVID Moonshot](https://covid.postera.ai/covid) project, as well as the initial experimental data created as part of the project.

The data has been split into many different folders and files, in order to ease the organization of such a large number of designs made through several CROs.

There exists one "master file" containing most of the non-experimental information: [covid_submissions_all_info.csv](https://github.com/mc-robinson/COVID_moonshot_submissions/blob/master/covid_submissions_all_info.csv) This file is described in much greater detail at the end of this README.

There also exist several folders:
- [data_for_CDD/](https://github.com/mc-robinson/COVID_moonshot_submissions/tree/master/data_for_CDD): prepares the data in a format suitable for upload to the CDD (Collaborative Drug Discovery) vault
- [orders_data/](https://github.com/mc-robinson/COVID_moonshot_submissions): The compounds that have been ordered up to this point.
- [shipments_data/](https://github.com/mc-robinson/COVID_moonshot_submissions/tree/master/shipments_data) the compounds that have been recevied up to this point.

Some more detail on [covid_submissions_all_info.csv](https://github.com/mc-robinson/COVID_moonshot_submissions/blob/master/covid_submissions_all_info.csv), and its many fields.

This is the main file containing all of the submitted molecules. It contains quite a few fields described below:
- **SMILES**: the SMILES string as canonicalized by RDKit and standardized by https://github.com/chembl/ChEMBL_Structure_Pipeline.
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
- **Covalent Fragment**: if one of the fragments cited as inpiration is a covalent fragment
- **covalent_warhead**: if *at least one* covalent warhead moiety is present in the designed molecule.
- **Acrylamide**: If the acrylamide moiety is in the molecule
- **Acrylamide Adduct**: If the acrylamide adduct moiety is in the molecule
- **Chloroacetamide**: If the chloroacetamide moiety is in the molecule
- **Chloroacetamide Adduct**: If the chloroacetamide adduct moiety is in the molecule 
- **Vinylsulfonamide**: If the vinylsulfonamide moiety is in the molecule
- **Vinylsulfonamide Adduct**: If the vinylsulfonamide adduct moiety is in the molecule
- **Nitrile**: If the nitrle moiety is in the molecule
- **Nitrile Adduct**: If the nitrle adduct moiety is in the molecule
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

## Assay Details ##

### Fluorescence MPro assay: ###

Compounds were seeded into assay-ready plates (Greiner 384 low volume 784900) using an Echo 555 acoustic dispenser, and DMSO was back-filled for a uniform concentration in assay plates (maximum 1%). Screening assays were performed in duplicate at 20 µM and 50 µM. Hits of greater than 50% inhibition at 50 µM were confirmed by dose response assays. Reagents for Mpro assay reagents were dispensed into the assay plate in 10 µl volumes for a final of 20 µl. Final reaction concentrations were 20 mM HEPES pH=7.3, 1mM TCEP, 50 mM NaCl, 0.01% Tween-20, 10% glycerol, 5 nM Mpro, 375 nM fluorogenic peptide substrate ([5-FAM]-AVLQSGFR-[Lys(Dabcyl)]-K-amide). Mpro was pre-incubated for 15 minutes at room temperature with compound before addition of substrate. Protease reaction was measured continuously in a BMG Pherastar FS with a 480/520 ex/em filter set.

### RapidFire MPro assay: ###

Inhibitor compounds are dispensed into 384-well plates using the ECHO dispenser (DMSO concentration < 1%). Enzyme solution, containing 20 nM Mpro, 20 mM HEPES, pH 7.5 and 50 mM NaCl, is added to each well and incubated with the inhibitor for 15 min at RT. The reaction is initiated with the addition of 2.0 μM substrate (TSAVLQSGFRK, custom synthesized in Schofield group). After 10 min the reaction is quenched with 10% formic acid and injected into an Agilent RapidFire LC-MS system. Data analysis are done with PRISM and CDD.  All compounds are triaged by testing calculating the % inhibition at 5 and 50 μM final concentration.  Dose response curves are done with 11 datapoints range of 100 – 0.0017 μM inhibitor.

### STD NMR ### 

*Instruments*
- Oxford Biochemistry 950 MHz NMR (Oxford Instruments magnet, Bruker AvanceIII console, 5mm TCI cryoprobe with tuning accessory, SampleJet temperature-controlled sample changer)
- Eppendorf 5810R centrifuge with deepwell plate rotorEppendorf 5355 thermomixer
- Eppendorf Repeater Plus repeating pipetteGilson pipettes

*Sample preparation procedure*
Final concentrations in NMR samples:
- 5% v/v D6-DMSO
- 10% v/v D2O 
- 10 uM Mpro
- 90 uM ligand

- Centrifuge ligand-containing 96-well plates at 3 500 rpm, 4 oC, for 5 minutes.
- Dissolve the ligand adding 7.75 ul of 100% v/v D6-DMSO to each well of the plate. Dispense 147.25 ul of PBS supplemented with D2O (ligand-only control experiments) or PBS supplemented with D2O and Mpro (STD experiments).
- Seal 96-well plate with adhesive aluminium foils and shake at 500 rpm at RT for 30 min.
- Centrifuge 96-well plates at 3 500 rpm, 4 oC, for 5 minutes.
- Transfer (in multiple, small steps) at least 140 μl of each ligand-containing solution to racked, 3 mm NMR tubes, using gel-loading tips.
- Seal NMR tubes with pom-balls, record position and barcode per sample

*NMR procedure*
- Pre-set temperature of NMR instrument and sample changer (SC) to 10 oC. Equilibrate for at least 30 min.
- Insert NMR tube racks to SC and record the rack position. Set the SC to 3 mm shuttle mode.
- Using TopSpin (currently, version 3.6.1, topspin commands or pulse sequences in parenthesis.
- Insert first sample in magnet with SC (sx). Let equilibrate for at least 4 min.
- Tune/Match ATM probe manually (atmm).Read previous shim-set, optimised for 3 mm tubes (rsh)
- Lock on D2O frequency using 90% H2O / 10% D2O preset (lock).
- Perform 3D gradient shim followed by temperature-optimised gradient shim (topshim 3d; tso).Optimise lock phase (autophase).
- Calibrate 90o 1H pulse length, record H2O frequency (zg pulse sequence).
- Enter new 1H pulse length and H2O frequency in otherwise optimised STD-NMR parameter set for Mpro (stddiffesgp.3 pulse sequence).
- Save experiment parameters (wpar) and sample shims (wsh).
- Using IconNMR Automation (currently, version 5.0.9 build 47):
- Setup automation experiment using saved experiment parameters and sample shims.
- For each new sample allow 4 min of temperature equilibration.
- For each new sample perform automatic tuning/matching (atma), locking using 90% H2O / 10% D2O preset (lock), gradient shimming (topshim), receiver gain adjustment (rga).
- Each new sample takes ~1 h total time, including overheads.

*Data analysis*
Using TopSpin (currently, version 3.6.1, topspin commands in parenthesis):
- Deconvolute the STD and reference spectra (stdsplit) into 2 directories. Reference frequency is #2.
- In the reference spectrum, apply window function, FT and phase the spectra (efp; apk). Record the phase correction values PHC0 and PHC1.
- Copy PHC0 and PHC1 values to STD spectrum, apply window function, FT and phase (efp).
- Manually compare the reference spectrum to presumed ligand chemistry. Note impurities and/or concerns. 
- Integrate the areas of DMSO and ligand peaks in the reference spectrum. Derive ratio of integrated areas between DMSO and single 1H of ligand as metric of relative ligand concentration.
- Integrate the areas of ligand peaks in STD spectrum. Derive the ratio (multiplied by 103) of ligand peak areas between STD and reference spectrum as metric of ligand binding to Mpro.

### High-throughput solubility threshold measurement ###

Nephelometry-based solubility assay to define threshold compound solubility at the concentrations of 20 µM and 100 µM in PBS solution. The compound solubilities are normalized to Deoxyfluorouridine (100% soluble) and ondansetron (0% solubility).

To put the numbers in perspective, the expected relative solubility values are:
- High >0.8
- Mid 0.6-0.8
- Low <0.6

To perform the screening, 245 µL aliquots of PBS buffer were added to each well of 96-well
microplates with clear flat bottom. Plates with the buffer were subjected to optical integrity
inspection using Nephelostar. Total scanning time for one 96-well plate was 3 min.

The pass criterion was set as the background signal in any of the scanned wells below 25 RNU, thus making optical quality of the plates satisfactory for the assay. In case of an excessively high background signal in any wells of the test-plates, those plates/wells were excluded from the study (data not shown).

Compounds were obtained from Enamine repository as solids formatted in polypropylene, round
bottom blank tubes, in latch boxes (Matrix #4271). Compounds were dissolved in DMSO at
50 mM, incubated at 24-26°C for 8 hours, shaken for 1 hour at 1800 rpm using high-speed
microplate shaker Illumina, then incubated at 24-26°C for 14 hours. Intermediate DMSO solutions
of the tested compounds were prepared to 2 mM (for the tested concentration 20 µM) and 8 mM
(for the tested concentration 100 µM). Ondansetron and DOFU were added to the final plates in
100 mM concentration to get 1 mM and 2 mM final concentrations in the assay plates.

To prepare the test plates, 2.5µL aliquots of DMSO solutions of the tested compounds, reference
compounds or pure DMSO were transferred from the polypropylene tubes to the corresponding
wells of 96-well plates with PBS buffer using Plate Mate, according to the plate map (Fig. 1). Thus, final volumes in the test plates were brought to 247.5 µL, resulting in concentrations of the compounds in test wells of 20 µM and 1% DMSO, correspondingly. Then, turbidity of the solutions was immediately scanned for each well. 
















