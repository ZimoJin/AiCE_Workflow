from pathlib import Path

input_dir = Path("input")
output_dir = Path("output/single_result")
output_dir.mkdir(parents=True, exist_ok=True)

ss_dict = {"H", "B", "E", "G", "I", "T", "S", "P"}  # structured codes

dssp_files = list(input_dir.glob("*.dssp"))
if not dssp_files:
    print("❌ No .dssp files found in input/")
    exit(1)

for dssp_file in dssp_files:
    prefix = dssp_file.stem
    output_file = output_dir / f"{prefix}.ss"

    lines = dssp_file.read_text().splitlines()
    try:
        start_idx = next(i for i, l in enumerate(lines) if l.startswith("  #"))
    except StopIteration:
        print(f"⚠️  Skipping {dssp_file.name} — no header line found.")
        continue

    records = lines[start_idx+1:]

    with output_file.open("w") as out:
        for line in records:
            try:
                res_num = int(line[5:10].strip())
                aa = line[13].strip()
                ss = line[16].strip()
                if ss not in ss_dict:
                    ss = "C"
                out.write(f"{res_num}\t{aa}\t{ss}\n")
            except Exception:
                continue

    print(f"✅ Converted {dssp_file.name} ➝ {output_file}")
