import gzip

# Load your 8,000 whitelist lines into a set for fast lookup
print("Loading whitelist...")
with open("barcode_list.txt", "r") as f:
    whitelist = set(line.strip() for line in f)

r1_in = "./fastq/ScHiC1_R1.fastq.gz"
r2_in = "./fastq/ScHiC1_R2.fastq.gz"
r1_out = "./fastq/trimmed_1bc_ScHiC1_R1.fastq.gz"
r2_out = "./fastq/trimmed_1bc_ScHiC1_R2.fastq.gz"

total_reads = 0
recovered_reads = 0

print("Processing FASTQ files via sliding-window barcode matching...")
with gzip.open(r1_in, "rt") as f1, gzip.open(r2_in, "rt") as f2, \
     gzip.open(r1_out, "wt") as out1, gzip.open(r2_out, "wt") as out2:
    
    while True:
        h1, s1, b1, q1 = f1.readline(), f1.readline(), f1.readline(), f1.readline()
        h2, s2, b2, q2 = f2.readline(), f2.readline(), f2.readline(), f2.readline()
        if not h1: break
        
        total_reads += 1
        s1, s2 = s1.strip(), s2.strip()
        h1, h2, q1, q2 = h1.strip(), h2.strip(), q1.strip(), q2.strip()
        
        extracted_bc = None
        clean_s1, clean_s2 = s1, s2
        clean_q1, clean_q2 = q1, q2
        
        # Look for a 30bp sequence window anywhere inside Read 1 or Read 2 
        # that matches a sequence in your whitelist
        for i in range(len(s1) - 30):
            window = s1[i:i+30]
            if window in whitelist:
                extracted_bc = window
                # Trim the barcode completely out of the genomic sequence
                clean_s1 = s1[:i] + s1[i+30:]
                clean_q1 = q1[:i] + q1[i+30:]
                break
                
        if not extracted_bc:
            for i in range(len(s2) - 30):
                window = s2[i:i+30]
                if window in whitelist:
                    extracted_bc = window
                    clean_s2 = s2[:i] + s2[i+30:]
                    clean_q2 = q2[:i] + q2[i+30:]
                    break

        # If a barcode was found (even if truncated downstream), save the pair
        if extracted_bc:
            recovered_reads += 1
            new_h1 = f"{h1.split()[0]} CB:Z:{extracted_bc}"
            new_h2 = f"{h2.split()[0]} CB:Z:{extracted_bc}"
            out1.write(f"{new_h1}\n{clean_s1}\n+\n{clean_q1}\n")
            out2.write(f"{new_h2}\n{clean_s2}\n+\n{clean_q2}\n")

print(f"\n--- Barcode Recovery Summary ---")
print(f"Total read pairs evaluated: {total_reads:,}")
print(f"Read pairs with found barcodes: {recovered_reads:,} ({(recovered_reads/total_reads)*100:.2f}%)")
