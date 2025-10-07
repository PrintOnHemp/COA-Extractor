COA EXTRACTOR
Small utility to extract fields from lab COA PDFs and normalize them into a CSV.


VERSION
- 0.2.7
- 2025-10-07 19:40


USAGE
Run the Command Line Interface to process the `input/` folder and write to `output/`:
```bash
python3 coa_extractor.py
```


REQUIREMENTS
See `requirements.txt` for the Python dependencies.


CREDITS
- Vibe Code Contributions by Matt Glyer and Austin Smith using ChatGPT5.


VERSION HISTORY
- 0.2.7 (2025-10-07 19:40): Implemented global 'round half up' rounding rule to replace Python's default banker’s rounding, ensuring consistent human-expected rounding for all numeric and percentage outputs.
- 0.2.6 (2025-10-07 17:30): Added global date normalization rule; all extracted dates now formatted as mm/dd/yyyy for consistency across labs.
- 0.2.5 (2025-10-07 16:55): Strengthened Green Leaf Lab terpene parsing by splitting on the full LOD–LOQ–Units pattern (`0.002 0.005 mg/g`) to correctly capture whole-number results such as β-Caryophyllene 5.
- 0.2.4 (2025-10-07 16:25): Hardened Green Leaf Lab terpene parsing to accept integer results (e.g., “5”) as valid mg/g values, including cases where `mg/g` splits across a newline. Improved analyte matching to include Greek letters and normalized whitespace.
- 0.2.3 (2025-10-07 16:08): Fixed terpene parsing for whole-number results (e.g., β-Caryophyllene = 5) in Green Leaf Lab data.
- 0.2.2 (2025-10-07 15:42): Updated formatting for all percentages to display two decimal places and confirmed accurate ChemHistory parsing improvements.