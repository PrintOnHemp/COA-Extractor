COA EXTRACTOR
Small utility to extract fields from lab COA PDFs and normalize them into a CSV.


VERSION
- 0.2.14
- 2025-10-08 17:20


USAGE
Run the Command Line Interface to process the `input/` folder and write to `output/`:
```bash
python3 coa_extractor.py
```


REQUIREMENTS
See `requirements.txt` for the Python dependencies.


CREDITS
- Vibe Code Contributions by Matt Glyer and Austin Smith using ChatGPT5.


- 0.2.14 (2025-10-08 17:20): Reordered CSV output columns to match the user-defined final structure: METRC Batch, Farm Name, Farm License, Strain Name, Testing Lab, Tested Date, Harvested Date, Total THC %, Total CBD %, Total Terpenes %, Total Cannabinoids %, Top Terpene 1–5 (and %), Source File.
- 0.2.12 (2025-10-08 17:10): Updated CSV output column order to match latest COA Parsed Data Template (2025_10_08).
- 0.2.11 (2025-10-08 16:00): Refined ChemHistory Tested Date extraction; removed fallback that used 'Date Accepted' and retained only relevant regex patterns targeting dates near 'Cannabinoids Pass'.
- 0.2.10 (2025-10-08 15:25): Enhanced Green Leaf Lab Tested Date extraction to handle different dash types (–, —) and newline-separated date formats following 'Lab Director'.
- 0.2.9 (2025-10-08 14:45): Removed fallback Tested Date extraction for Green Leaf Lab; line will remain empty if ‘Lab Director -’ date is not found.
- 0.2.8 (2025-10-08 14:10): Updated Green Leaf Lab Tested Date extraction to pull from the Lab Director signature line instead of 'Date Accepted' for improved accuracy.
- 0.2.7 (2025-10-07 19:40): Implemented global 'round half up' rounding rule to replace Python's default banker’s rounding, ensuring consistent human-expected rounding for all numeric and percentage outputs.
- 0.2.6 (2025-10-07 17:30): Added global date normalization rule; all extracted dates now formatted as mm/dd/yyyy for consistency across labs.
- 0.2.5 (2025-10-07 16:55): Strengthened Green Leaf Lab terpene parsing by splitting on the full LOD–LOQ–Units pattern (`0.002 0.005 mg/g`) to correctly capture whole-number results such as β-Caryophyllene 5.
- 0.2.4 (2025-10-07 16:25): Hardened Green Leaf Lab terpene parsing to accept integer results (e.g., “5”) as valid mg/g values, including cases where `mg/g` splits across a newline. Improved analyte matching to include Greek letters and normalized whitespace.
- 0.2.3 (2025-10-07 16:08): Fixed terpene parsing for whole-number results (e.g., β-Caryophyllene = 5) in Green Leaf Lab data.
- 0.2.2 (2025-10-07 15:42): Updated formatting for all percentages to display two decimal places and confirmed accurate ChemHistory parsing improvements.