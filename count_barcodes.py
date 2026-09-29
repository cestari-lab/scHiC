import gzip
import collections

# Input file path - Change this to match your actual file name
fastq_file = "/mnt/h/Analysis_2/scHi-C/fastq/trimmed_1bc_ScHiC1_R1.fastq.gz"  # or "chromap_ready_R1.fastq.gz"

# Dictionary to hold barcode counts
barcode_counts = collections.Counter()

print("Scanning FASTQ file and talling barcodes...")

# Open file automatically handling gzip if compressed
open_func = gzip.open if fastq_file.endswith(".gz") else open
mode = "rt" if fastq_file.endswith(".gz") else "r"

total_reads = 0

with open_func(fastq_file, mode) as f:
    while True:
        # FASTQ files are structured in 4-line blocks
        header = f.readline()
        seq = f.readline()
        plus = f.readline()
        qual = f.readline()
        
        if not header:
            break  # End of file reached
            
        total_reads += 1
        
        # Look for the cell barcode identifier tag
        if "CB:Z:" in header:
            # Extract everything after the tag, strip whitespace or trailing elements
            barcode = header.strip().split("CB:Z:")[-1]
            barcode_counts[barcode] += 1
        else:
            barcode_counts["NO_BARCODE_TAG"] += 1

# Convert counts to a sorted list (Highest abundance first)
sorted_groups = barcode_counts.most_common()

# --- Print the Final Summary Report ---
print("\n" + "="*50)
print("             SINGLE-CELL BARCODE REPORT            ")
print("="*50)
print(f"Total sequences processed : {total_reads:,}")
print(f"Unique cell groups found  : {len(barcode_counts):,}")
print("-"*50)
print(f"{'Rank':<6}{'Cell Barcode Sequence':<34}{'Read Count':<12}{'Percentage':<10}")
print("-"*50)

for rank, (barcode, count) in enumerate(sorted_groups, start=1):
    percentage = (count / total_reads) * 100
    print(f"{rank:<6}{barcode:<34}{count:<12,}{percentage:<10.2f}%")
    
    # Optional: Stop printing to terminal screen after the top 20 cells to save space
    if rank >= 20:
        print(f"... and {len(barcode_counts) - 20:,} more cell groups.")
        break
        
print("="*50)

# --- Save Full Detailed Table to a CSV File ---
csv_output = "/mnt/h/Analysis_2/scHi-C/fastq/barcode_summary_report.csv"
with open(csv_output, "w") as out:
    out.write("Rank,Barcode,Read_Count,Percentage\n")
    for rank, (barcode, count) in enumerate(sorted_groups, start=1):
        percentage = (count / total_reads) * 100
        out.write(f"{rank},{barcode},{count},{percentage:.4f}\n")

print(f"\nFull detailed summary matrix saved to: {csv_output}")
