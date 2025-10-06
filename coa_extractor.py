import os, re, pdfplumber, pandas as pd
from datetime import datetime

def read_first_n_pages_text(pdf_path, n=3):
    texts = []
    with pdfplumber.open(pdf_path) as pdf:
        for i, page in enumerate(pdf.pages[:n]):
            texts.append(page.extract_text() or "")
    return "\n".join(texts)

def clean_num(val):
    if val is None: return None
    s = str(val).strip().replace("%","").replace("ppm","").replace("mg/g","")
    if "< LOQ" in s or "<LOQ" in s:
        return "<LOQ"
    if s.upper() in {"ND", "N.D."}:
        return "ND"
    import re
    m = re.search(r"[-+]?\d*\.?\d+", s)
    return round(float(m.group(0)), 2) if m else None

def extract_top_terpenes(text):
    import re
    terpenes = []
    # Limit the text to the terpene analysis section
    m_section = re.search(r"Terpene Analysis by GCMS(.*?)Total Terpenes", text, re.S | re.I)
    if not m_section:
        return [(None, None), (None, None), (None, None), (None, None), (None, None)]
    terp_text = m_section.group(1)
    lines = terp_text.strip().splitlines()
    # Remove the first two lines if they contain date/time extraction and analysis info
    if len(lines) > 2:
        lines = lines[2:]
    terp_text = "\n".join(lines)
    # Split text into lines and then by mg/g to isolate each terpene entry
    for line in terp_text.splitlines():
        # Break line into individual terpene chunks at 'mg/g'
        chunks = re.split(r"\s*mg/g\s*", line)
        for chunk in chunks:
            if not chunk.strip():
                continue
            m = re.match(r"([A-Za-z0-9\-\(\)\s]+)\s+([0-9\.<>NDnd]+)", chunk.strip())
            if not m:
                continue
            analyte = m.group(1).strip()
            greek_map = {
                r"\balpha\b": "α",
                r"\bbeta\b": "β",
                r"\bgamma\b": "γ",
                r"\bdelta\b": "δ",
                r"\bepsilon\b": "ε",
                r"\bzeta\b": "ζ",
                r"\beta\b": "β",
                r"\btheta\b": "θ",
                r"\blambda\b": "λ",
                r"\bmu\b": "μ",
                r"\bnu\b": "ν",
                r"\bpi\b": "π",
                r"\brho\b": "ρ",
                r"\bsigma\b": "σ",
                r"\btau\b": "τ",
                r"\bphi\b": "φ",
                r"\bchi\b": "χ",
                r"\bpsi\b": "ψ",
                r"\bomega\b": "ω"
            }
            for k, v in greek_map.items():
                analyte = re.sub(k, v, analyte, flags=re.I)
            result_str = m.group(2).strip()
            if "nd" in result_str.lower() or "total" in analyte.lower():
                continue
            if "<" in result_str:
                terpenes.append((analyte, "<LOQ"))
                continue
            try:
                # Convert mg/g to percentage for Green Leaf Lab by dividing by 10
                result = round(float(result_str) / 10, 2)
                terpenes.append((analyte, result))
            except:
                continue
    numeric_terpenes = [t for t in terpenes if isinstance(t[1], (int, float))]
    loq_terpenes = [t for t in terpenes if t[1] == "<LOQ"]
    numeric_terpenes.sort(key=lambda x: x[1], reverse=True)
    terpenes = numeric_terpenes + loq_terpenes
    top5 = []
    for name, val in terpenes[:5]:
        if isinstance(val, (int, float)):
            top5.append((name, f"{val}%"))
        elif val == "<LOQ":
            top5.append((name, "<LOQ"))
        else:
            top5.append((name, None))
    while len(top5) < 5:
        top5.append((None, None))
    return top5

def extract_top_terpenes_chemhistory(text):
    import re
    terpenes = []
    m_section = re.search(r"Analyte Mass Mass LOQ.*?% mg/g % % mg/g %(.*?)(?:Primary Aromas)", text, re.S|re.I)
    if not m_section:
        return [(None, None)] * 5
    terp_text = m_section.group(1).strip()
    lines = terp_text.splitlines()
    for line in lines:
        # Each analyte in ChemHistory has: Name  <val/% or <LOQ>  <mg/g or <LOQ>  <LOQ-limit-number>
        # The trailing LOQ-limit number (e.g., 0.02) was being left behind and then prefixed to the next analyte.
        # Fix: (1) require analyte to START with a letter (Latin or Greek) so numbers can't prefix names,
        #      (2) consume the trailing LOQ-limit token, if present.
        matches = re.findall(
            r"([A-Za-z\(\)\-\u03B1-\u03C9\u0391-\u03A9][A-Za-z0-9\-\(\)\s\u03B1-\u03C9\u0391-\u03A9]*?)"
            r"\s+((?:<\s*LOQ)|(?:[0-9]+\.[0-9]+))"
            r"\s+((?:<\s*LOQ)|(?:[0-9]+\.[0-9]+))"
            r"\s+(?:<\s*LOQ|[0-9]+\.[0-9]+)?",
            line
        )
        for analyte, first_val, second_val in matches:
            analyte = analyte.strip()
            if "<" in first_val.upper():
                terpenes.append((analyte, "<LOQ"))
                continue
            try:
                value = round(float(first_val), 2)
                terpenes.append((analyte, value))
            except:
                continue
    numeric_terpenes = [t for t in terpenes if isinstance(t[1], (int, float))]
    loq_terpenes = [t for t in terpenes if t[1] == "<LOQ"]
    numeric_terpenes.sort(key=lambda x: x[1], reverse=True)
    terpenes = numeric_terpenes + loq_terpenes
    top5 = []
    for name, val in terpenes[:5]:
        if isinstance(val, (int, float)):
            top5.append((name, f"{val}%"))
        elif val == "<LOQ":
            top5.append((name, "<LOQ"))
        else:
            top5.append((name, None))
    while len(top5) < 5:
        top5.append((None, None))
    return top5

def extract_confident_lims(text):
    out = {"Testing Lab":"ChemHistory"}
    import re
    m = re.search(r"Sample:\s*([A-Z0-9\.\-]+)", text, re.I)
    out["Strain Name"] = None
    m_strain = re.search(r"Strain:\s*([^\n]+)", text, re.I)
    if m_strain:
        out["Strain Name"] = m_strain.group(1).strip()
    else:
        # fallback: try to find strain below Sample:
        sample_block = None
        sample_match = re.search(r"Sample:\s*([A-Z0-9\.\-]+)(.*?)(?:\n\n|\Z)", text, re.I | re.S)
        if sample_match:
            sample_block = sample_match.group(2)
            m_strain_fallback = re.search(r"([^\n]+)", sample_block.strip())
            if m_strain_fallback:
                out["Strain Name"] = m_strain_fallback.group(1).strip()
    # Extract Farm Name as the text between "1 of N" and "Sample:"
    m_farm = re.search(r"1 of \d+\s*\n([^\n]+)\s*Sample:", text, re.I)
    out["Farm Name"] = m_farm.group(1).strip() if m_farm else None
    # Broaden license extraction to letters + numbers
    m_license = re.search(r"Lic[^0-9A-Z]*([0-9A-Z\-]+)", text, re.I)
    if not m_license:
        m_license = re.search(r"Lic\s*#\s*([0-9A-Z\-]+)", text, re.I)
    out["Farm License"] = m_license.group(1).strip() if m_license else None

    # Extract cannabinoid percentages after 3rd PASS
    pass_matches = [m.start() for m in re.finditer(r"PASS", text, re.I)]
    out["Total THC %"], out["Total CBD %"], out["Total Cannabinoids %"] = None, None, None
    if len(pass_matches) >= 3:
        segment = text[pass_matches[2]:]
        percents = re.findall(r"([0-9]+\.[0-9]+)\s*%", segment)
        if len(percents) >= 1:
            out["Total THC %"] = f"{round(float(percents[0]), 2)}%"
        if len(percents) >= 2:
            out["Total CBD %"] = f"{round(float(percents[1]), 2)}%"
        if len(percents) >= 3:
            out["Total Cannabinoids %"] = f"{round(float(percents[2]), 2)}%"

    # Only capture the number if it appears on the same line as "Moisture Activity"
    terp_match = re.search(r"Moisture Activity[^\n]*?([0-9]+\.[0-9]+)%?", text, re.I)
    out["Total Terpenes %"] = f"{round(float(terp_match.group(1)), 2)}%" if terp_match else None

    m = re.search(r"Harvest/Production Date:\s*([0-9/]{2}/[0-9/]{2}/[0-9]{4})", text)
    if not m:
        m = re.search(r"Harvest/Prod\. Date:\s*([0-9\./\-]+)", text)
    out["Harvested Date"] = m.group(1) if m else None

    # Extract METRC Batch
    m_metrc = re.search(r"METRC Batch:\s*([A-Z0-9]+)", text, re.I)
    out["METRC Batch"] = m_metrc.group(1) if m_metrc else None

    # Improved Tested Date extraction: look for Cannabinoids section date first
    m_tested = re.search(r"Cannabinoids.*?\n\s*([0-9]{2}/[0-9]{2}/[0-9]{4})", text, re.I)
    if not m_tested:
        # fallback to any standalone date after "Cannabinoids"
        m_tested = re.search(r"Cannabinoids.*?(?:Pass)?\s*\n.*?([0-9]{2}/[0-9]{2}/[0-9]{4})", text, re.I|re.S)
    if not m_tested:
        m_tested = re.search(r"Date Accepted:\s*([0-9/]{2}/[0-9/]{2}/[0-9]{2,4})", text)
    out["Tested Date"] = m_tested.group(1) if m_tested else None
    out["Source File"] = None

    top5 = extract_top_terpenes_chemhistory(text)
    print("ChemHistory Top Terpenes:", top5)
    out["Top Terpene 1"], out["Top Terpene 1 %"] = top5[0]
    out["Top Terpene 2"], out["Top Terpene 2 %"] = top5[1]
    out["Top Terpene 3"], out["Top Terpene 3 %"] = top5[2]
    out["Top Terpene 4"], out["Top Terpene 4 %"] = top5[3]
    out["Top Terpene 5"], out["Top Terpene 5 %"] = top5[4]

    # No terpene extraction for ChemHistory/Confident LIMS (Top Terpene 1/2/3 not set)

    return out

def extract_greenleaf(text):
    import re
    out = {"Testing Lab":"Green Leaf Lab"}
    m = re.search(r"Sample ID:\s*([A-Z0-9\-]+)", text, re.I)
    out["Strain Name"] = None
    m_strain = re.search(r"(?m)^(.+?)\s*\nSample ID:", text)
    if m_strain:
        out["Strain Name"] = m_strain.group(1).strip()
    else:
        m_strain = re.search(r"\n([A-Za-z0-9#\-\s]+)\nMatrix:\s*Useable Marijuana", text)
        if m_strain:
            out["Strain Name"] = m_strain.group(1).strip()
    m_client = re.search(r"Client:\s*([^\n]+)", text, re.I)
    if not m_client:
        m_client = re.search(r"Producer:\s*([^\n]+)", text, re.I)
    out["Farm Name"] = m_client.group(1).strip() if m_client else None
    if not out["Farm Name"]:
        m_farm = re.search(r"Harvest/Prod\. Date:\s*[0-9\./-]+\s+([A-Za-z0-9 ,.&]+)", text)
        if m_farm:
            farm_name = m_farm.group(1).strip()
            farm_name = re.sub(r"\s+\d.*$", "", farm_name)  # remove trailing numbers/codes
            out["Farm Name"] = farm_name
    out["Farm License"] = None

    out["Total THC %"] = None
    m_loq = re.search(r"Total THC\s*:\s*<\s*LOQ", text, re.I)
    if m_loq:
        out["Total THC %"] = "<LOQ"
    else:
        m = re.search(r"Total THC\s*:\s*([0-9\.]+)\s*%", text, re.I)
        if m:
            out["Total THC %"] = f"{round(float(m.group(1)), 2)}%"

    out["Total CBD %"] = None
    m_loq = re.search(r"Total CBD\s*:\s*<\s*LOQ", text, re.I)
    if m_loq:
        out["Total CBD %"] = "<LOQ"
    else:
        m = re.search(r"Total CBD\s*:\s*([0-9\.]+)\s*%", text, re.I)
        if m:
            out["Total CBD %"] = f"{round(float(m.group(1)), 2)}%"

    out["Total Terpenes %"] = None
    m_loq = re.search(r"Total Terpenes\s*:\s*<\s*LOQ", text, re.I)
    if m_loq:
        out["Total Terpenes %"] = "<LOQ"
    else:
        m = re.search(r"Total Terpenes\s*:\s*([0-9\.]+)\s*%", text, re.I)
        if m:
            out["Total Terpenes %"] = f"{round(float(m.group(1)), 2)}%"

    # Extract Total Cannabinoids % as the first number after "Total Cannabinoids" (percent sign may not be present)
    out["Total Cannabinoids %"] = None
    m = re.search(r"Total Cannabinoids\s+([0-9]+\.[0-9]+)", text)
    if m:
        out["Total Cannabinoids %"] = f"{round(float(m.group(1)), 2)}%"
    m = re.search(r"Harvest/Prod\. Date:\s*([0-9\./]+)", text)
    if m:
        date_val = m.group(1).strip()
        # Convert dots to slashes for consistency
        date_val = date_val.replace(".", "/")
        out["Harvested Date"] = date_val
    else:
        out["Harvested Date"] = None
    m = re.search(r"Date Accepted:\s*([0-9/]{2}/[0-9/]{2}/[0-9]{2,4})", text)
    out["Tested Date"] = m.group(1) if m else None
    # Extract METRC Batch or Source ID
    m_metrc = re.search(r"(?:METRC Batch|Source ID):\s*([A-Z0-9]+)", text, re.I)
    out["METRC Batch"] = m_metrc.group(1) if m_metrc else None
    out["Source File"] = None

    top5 = extract_top_terpenes(text)
    out["Top Terpene 1"], out["Top Terpene 1 %"] = top5[0]
    out["Top Terpene 2"], out["Top Terpene 2 %"] = top5[1]
    out["Top Terpene 3"], out["Top Terpene 3 %"] = top5[2]
    out["Top Terpene 4"], out["Top Terpene 4 %"] = top5[3]
    out["Top Terpene 5"], out["Top Terpene 5 %"] = top5[4]

    return out

def detect_lab(text):
    tl = text.lower()
    if "green leaf lab" in tl: return "greenleaf"
    if "confident lims" in tl or "chemhistory" in tl: return "confident"
    return "unknown"

def extract_fields(pdf_path, pages=3):
    # First read enough to detect the lab
    text = read_first_n_pages_text(pdf_path, n=pages)
    lab = detect_lab(text)
    # If ChemHistory, read up to 5 pages for terpenes
    if lab == "confident":
        text = read_first_n_pages_text(pdf_path, n=5)
        out = extract_confident_lims(text)
    elif lab == "greenleaf":
        out = extract_greenleaf(text)
    # Removed generic extraction branch
    # else: out = extract_generic(text)
    # Set source file
    out["Source File"] = os.path.basename(pdf_path)
    return out

def process_folder(input_folder, output_csv="coa_unified.csv", pages=3):
    rows = []
    for name in os.listdir(input_folder):
        if not name.lower().endswith(".pdf"): continue
        path = os.path.join(input_folder, name)
        try:
            rows.append(extract_fields(path, pages=pages))
        except Exception as e:
            rows.append({"Source File": name, "Error": str(e)})
    cols = [
      "Strain Name","Farm Name","Farm License","Total Cannabinoids %","Total THC %","Total CBD %","Total Terpenes %",
      "Harvested Date","METRC Batch","Testing Lab","Tested Date","Source File",
      "Top Terpene 1","Top Terpene 1 %","Top Terpene 2","Top Terpene 2 %","Top Terpene 3","Top Terpene 3 %",
      "Top Terpene 4","Top Terpene 4 %","Top Terpene 5","Top Terpene 5 %"
    ]
    df = pd.DataFrame(rows, columns=cols)
    output_dir = os.path.dirname(os.path.abspath(output_csv))
    if output_dir and not os.path.isdir(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    df.to_csv(output_csv, index=False)
    return df

if __name__ == "__main__":
    import argparse
    script_dir = os.path.dirname(os.path.abspath(__file__))
    default_input = os.path.join(script_dir, "input")
    timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    default_output = os.path.join(script_dir, "output", f"coa_unified_{timestamp}.csv")
    p = argparse.ArgumentParser(description="COA PDF normalizer")
    p.add_argument(
        "input_folder",
        nargs="?",
        default=default_input,
        help=f"Folder containing lab PDFs (default: {default_input})",
    )
    p.add_argument("--pages", type=int, default=3, help="How many pages to parse per PDF")
    p.add_argument(
        "--out",
        default=default_output,
        help=f"Output CSV path (default: {default_output})",
    )
    args = p.parse_args()
    df = process_folder(args.input_folder, output_csv=args.out, pages=args.pages)
    print(df.head(10).to_string(index=False))
