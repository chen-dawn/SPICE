# Cloning Procedure

The primary library used for the screen was constructed using the **DC297** backbone (`DC297-Ef1a-3PieceGibsonVectorWEAK double BsmBI.gb`).

## Step 1: Library Insertion
The DC297 vector was digested with **BsmBI (Esp3I)**, and the library oligo pool was cloned into the plasmid using **Gibson assembly**. A representative library element is provided in `SAMPLE_library_element.gb`. The assembled library was transformed at sufficient coverage to ensure representation, generating the intermediate plasmid library `SAMPLE_final_plasmid_after_step1_cloning-DC297.gb`.

## Step 2: Constant Sequence Insertion
The constant sequence (`DCG197_ACTN1_Twist_BsmBI_dropin.gb`) was inserted into the pooled plasmid library using the **NEB Golden Gate BsmBI cloning protocol**. The resulting constructs were transformed into bacteria, yielding the final plasmid library `SAMPLE_final_plasmid_after_step2-cloning-DC297.gb`.
