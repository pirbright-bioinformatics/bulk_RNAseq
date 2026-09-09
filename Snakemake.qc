import os
import re
from pathlib import Path

configfile: "config.yaml"

# --- Robust Sample Discovery ---
def discover_samples(sampledir_path):
    sampledir = Path(sampledir_path)
    samples = {}
    
    # Common extensions for FASTQ files
    extensions = {".fq.gz", ".fastq.gz", ".fq", ".fastq"}
    
    # Patterns to identify paired-end mates (R1/1 and R2/2)
    pe_patterns = [
        r'^(.*)(_1)(\.fq|\.fastq)(\.gz)?$',
        r'^(.*)(_R1_001)(\.fq|\.fastq)(\.gz)?$',
        r'^(.*)(_R1)(\.fq|\.fastq)(\.gz)?$',
        r'^(.*)(\.1)(\.fq|\.fastq)(\.gz)?$',
    ]
    
    pe2_patterns = [
        r'^(.*)(_2)(\.fq|\.fastq)(\.gz)?$',
        r'^(.*)(_R2_001)(\.fq|\.fastq)(\.gz)?$',
        r'^(.*)(_R2)(\.fq|\.fastq)(\.gz)?$',
        r'^(.*)(\.2)(\.fq|\.fastq)(\.gz)?$',
    ]
    
    for root, dirs, files in os.walk(sampledir):
        for f in files:
            if not any(f.endswith(ext) for ext in extensions):
                continue
                
            is_pe = False
            prefix = f
            
            # Check if it matches R1 patterns
            for pattern in pe_patterns:
                m = re.match(pattern, f)
                if m:
                    prefix = m.group(1)
                    is_pe = True
                    break
            
            # If not matched as R1, check if it matches R2 patterns
            if not is_pe:
                for pattern in pe2_patterns:
                    m = re.match(pattern, f)
                    if m:
                        prefix = m.group(1)
                        is_pe = True
                        break
            
            # If still not matched, treat as single-end
            if not is_pe:
                for ext in sorted(extensions, key=len, reverse=True):
                    if f.endswith(ext):
                        prefix = f[:-len(ext)]
                        break
            
            # Calculate relative path to maintain directory structure
            rel_dir = Path(root).relative_to(sampledir)
            
            # De-duplicate folder name if it matches prefix (e.g., 4100/4100_1.fq.gz -> 4100)
            if rel_dir.name == prefix:
                sample_id = str(rel_dir)
            else:
                sample_id = str(rel_dir / prefix) if str(rel_dir) != "." else prefix
                
            if sample_id not in samples:
                samples[sample_id] = []
                
            abs_path = os.path.join(root, f)
            if abs_path not in samples[sample_id]:
                samples[sample_id].append(abs_path)
                
    # Sort files to ensure R1 comes before R2
    for sample_id in samples:
        samples[sample_id].sort()
        
    return samples

# Initialize parameters
SAMPLEDIR = Path(config["sampledir"])
SAMPLES_DICT = discover_samples(SAMPLEDIR)

PE_SAMPLES = [s for s, files in SAMPLES_DICT.items() if len(files) == 2]
SE_SAMPLES = [s for s, files in SAMPLES_DICT.items() if len(files) == 1]
ALL_SAMPLES = list(SAMPLES_DICT.keys())

print(f"Discovered {len(ALL_SAMPLES)} samples total.")
print(f"  Paired-End: {len(PE_SAMPLES)}")
print(f"  Single-End: {len(SE_SAMPLES)}")

MINLEN = config.get("trimmer", {}).get("minlen", 120)
OUTDIR = config["outdir"]
QCTYPE = config.get("qctype", "trimmed")
if QCTYPE == "": 
   QCTYPE = "trimmed"

rule all:
   input:
        # Raw FastQC
        fastqc_pe_1=expand(OUTDIR + "/fastqc/{ID}_1_fastqc.zip", ID=PE_SAMPLES),
        fastqc_pe_2=expand(OUTDIR + "/fastqc/{ID}_2_fastqc.zip", ID=PE_SAMPLES),
        fastqc_se=expand(OUTDIR + "/fastqc/{ID}_1_fastqc.zip", ID=SE_SAMPLES),
        
        # Trimmed FastQC
        fastqcp_pe_1=expand(OUTDIR + "/fastqc/{ID}/1."+QCTYPE+"_fastqc.zip", ID=PE_SAMPLES),
        fastqcp_pe_2=expand(OUTDIR + "/fastqc/{ID}/2."+QCTYPE+"_fastqc.zip", ID=PE_SAMPLES),
        fastqcp_se=expand(OUTDIR + "/fastqc/{ID}/1."+QCTYPE+"_fastqc.zip", ID=SE_SAMPLES),
        
        # Centralized Links
        OGlink_pe_1=expand(OUTDIR + "/fastqc/OG/{ID}/1_fastqc.zip", ID=PE_SAMPLES),
        OGlink_pe_2=expand(OUTDIR + "/fastqc/OG/{ID}/2_fastqc.zip", ID=PE_SAMPLES),
        OGlink_se=expand(OUTDIR + "/fastqc/OG/{ID}/1_fastqc.zip", ID=SE_SAMPLES),
        
        QCTYPElink_pe_1=expand(OUTDIR + "/fastqc/"+QCTYPE+"/{ID}/1_fastqc.zip", ID=PE_SAMPLES),
        QCTYPElink_pe_2=expand(OUTDIR + "/fastqc/"+QCTYPE+"/{ID}/2_fastqc.zip", ID=PE_SAMPLES),
        QCTYPElink_se=expand(OUTDIR + "/fastqc/"+QCTYPE+"/{ID}/1_fastqc.zip", ID=SE_SAMPLES),
   output:
        multiqc1=OUTDIR+"/fastqc/OG/QC_report_multiqc_report.html",
        multiqc2=OUTDIR+"/fastqc/"+QCTYPE +"/QC_report_multiqc_report.html"
   shell:
       """
       multiqc {OUTDIR}/fastqc/OG --interactive --title "QC_report" -f -o {OUTDIR}/fastqc/OG
       multiqc {OUTDIR}/fastqc/{QCTYPE} --interactive --title "QC_report" -f -o {OUTDIR}/fastqc/{QCTYPE}
       echo Done
       """

# --- PAIRED-END RULES ---

rule fastqc_pe:
   input:
        f1=lambda wildcards: SAMPLES_DICT[wildcards.ID][0],
        f2=lambda wildcards: SAMPLES_DICT[wildcards.ID][1],
   output:
        fastqc1=OUTDIR+"/fastqc/{ID}_1_fastqc.zip",
        fastqc2=OUTDIR+"/fastqc/{ID}_2_fastqc.zip"
   wildcard_constraints:
        ID="|".join(re.escape(s) for s in PE_SAMPLES) if PE_SAMPLES else "$^"
   shell:
       """
        mkdir -p $(dirname {output.fastqc1})
        
        tmp1=$(mktemp -d)
        fastqc {input.f1} -o $tmp1
        mv $tmp1/*_fastqc.zip {output.fastqc1}
        rm -rf $tmp1
        
        tmp2=$(mktemp -d)
        fastqc {input.f2} -o $tmp2
        mv $tmp2/*_fastqc.zip {output.fastqc2}
        rm -rf $tmp2
       """

rule fastqcp_pe:
   input:
        f1=OUTDIR +"/" +QCTYPE+"/{ID}_1.fq.gz",
        f2=OUTDIR +"/" +QCTYPE+"/{ID}_2.fq.gz"
   output:
        fastqc1=OUTDIR+"/fastqc/{ID}/1."+QCTYPE+"_fastqc.zip",
        fastqc2=OUTDIR+"/fastqc/{ID}/2."+QCTYPE+"_fastqc.zip"
   wildcard_constraints:
        ID="|".join(re.escape(s) for s in PE_SAMPLES) if PE_SAMPLES else "$^"
   shell:
       """
        mkdir -p $(dirname {output.fastqc1})
        
        tmp1=$(mktemp -d)
        fastqc {input.f1} -o $tmp1
        mv $tmp1/*_fastqc.zip {output.fastqc1}
        rm -rf $tmp1
        
        tmp2=$(mktemp -d)
        fastqc {input.f2} -o $tmp2
        mv $tmp2/*_fastqc.zip {output.fastqc2}
        rm -rf $tmp2
       """

rule fastp_pe:
   input:
        f1=lambda wildcards: SAMPLES_DICT[wildcards.ID][0],
        f2=lambda wildcards: SAMPLES_DICT[wildcards.ID][1],
   output:
        fastqc1=OUTDIR+"/"+QCTYPE+"/{ID}_1.fq",
        fastqc2=OUTDIR+"/"+QCTYPE+"/{ID}_2.fq",
   params:
        adapters=config.get("adapters", "adapters.fa"),
   log:
        stdout=OUTDIR +"/"+QCTYPE+ "/{ID}/log/stdout.txt",
        stderr=OUTDIR +"/"+QCTYPE+ "/{ID}/log/stderr.txt",
   threads: 8
   wildcard_constraints:
        ID="|".join(re.escape(s) for s in PE_SAMPLES) if PE_SAMPLES else "$^"
   shell:
       """
       mkdir -p $(dirname {output.fastqc1})/log

       if [ "{QCTYPE}" == "fastp" ];
       then
               if [ "{params.adapters}" != "" ] && [ -f "{params.adapters}" ];
               then
                       ADAPTER="--adapter_fasta {params.adapters}" 
               else
                       ADAPTER=""
               fi
               
               fastp --overrepresentation_analysis $ADAPTER --detect_adapter_for_pe --trim_poly_g \
                  --cut_mean_quality 20 --trim_poly_x --length_required 30 --thread {threads} \
                  --cut_window_size 4 --cut_front --cut_tail  \
                  --html $(dirname {output.fastqc1})/{QCTYPE}.html  \
                  --json $(dirname {output.fastqc1})/{QCTYPE}.json \
                  -i {input.f1}  -I {input.f2} -o {output.fastqc1}.tmp -O {output.fastqc2}.tmp > {log.stdout}  2>{log.stderr}

               cutadapt -a 'A{{10}}' -a 'T{{10}}' -A 'A{{10}}' -A 'T{{10}}' \
                  -o $(dirname {output.fastqc1})/fastp1c \
                  -p $(dirname {output.fastqc1})/fastp2c \
                  {output.fastqc1}.tmp {output.fastqc2}.tmp >> {log.stdout} 2>>{log.stderr}
                  
               rm {output.fastqc1}.tmp {output.fastqc2}.tmp
               
               /mnt/lustre/RDS-live/downing/FASTQC/fastq_quality_trimmer -Q 33 -t 30 -l {MINLEN} -i $(dirname {output.fastqc1})/fastp1c -o {output.fastqc1} -v >> {log.stdout} 2>>{log.stderr}
               /mnt/lustre/RDS-live/downing/FASTQC/fastq_quality_trimmer -Q 33 -t 30 -l {MINLEN} -i $(dirname {output.fastqc1})/fastp2c -o {output.fastqc2} -v >> {log.stdout} 2>>{log.stderr}
               
               rm $(dirname {output.fastqc1})/fastp1c $(dirname {output.fastqc1})/fastp2c
       else
               /mnt/lustre/RDS-live/downing/FASTQC/fastq_quality_trimmer -Q 33 -t 30 -l {MINLEN} -i <(zcat {input.f1}) -o {output.fastqc1} -v > {log.stdout} 2>{log.stderr}
               /mnt/lustre/RDS-live/downing/FASTQC/fastq_quality_trimmer -Q 33 -t 30 -l {MINLEN} -i <(zcat {input.f2}) -o {output.fastqc2} -v >> {log.stdout} 2>>{log.stderr}
       fi
       """

rule gzip_pe:
   input:
        fastqc1=OUTDIR + "/"+QCTYPE+"/{ID}_1.fq",
        fastqc2=OUTDIR + "/"+QCTYPE+"/{ID}_2.fq"
   output:
        gz1=OUTDIR + "/"+QCTYPE+"/{ID}_1.fq.gz",
        gz2=OUTDIR + "/"+QCTYPE+"/{ID}_2.fq.gz"
   wildcard_constraints:
        ID="|".join(re.escape(s) for s in PE_SAMPLES) if PE_SAMPLES else "$^"
   shell:
       """
       gzip {input.fastqc1}
       gzip {input.fastqc2}
       """

rule link_pe:
   input:
        fastqc1=OUTDIR+"/fastqc/{ID}_1_fastqc.zip",
        fastqc2=OUTDIR+"/fastqc/{ID}_2_fastqc.zip",
        f1=OUTDIR+"/fastqc/{ID}/1."+QCTYPE+"_fastqc.zip",
        f2=OUTDIR+"/fastqc/{ID}/2."+QCTYPE+"_fastqc.zip"
   output:
        OGln1=OUTDIR+"/fastqc/OG/{ID}/1_fastqc.zip",
        OGln2=OUTDIR+"/fastqc/OG/{ID}/2_fastqc.zip",
        FPln1=OUTDIR+"/fastqc/"+QCTYPE+"/{ID}/1_fastqc.zip",
        FPln2=OUTDIR+"/fastqc/"+QCTYPE+"/{ID}/2_fastqc.zip",
   wildcard_constraints:
        ID="|".join(re.escape(s) for s in PE_SAMPLES) if PE_SAMPLES else "$^"
   shell:
        """
        mkdir -p $(dirname {output.OGln1})
        mkdir -p $(dirname {output.FPln1})
        ln -sf $(readlink -f {input.fastqc1}) {output.OGln1}
        ln -sf $(readlink -f {input.fastqc2}) {output.OGln2}
        ln -sf $(readlink -f {input.f1}) {output.FPln1}
        ln -sf $(readlink -f {input.f2}) {output.FPln2}
        """

# --- SINGLE-END RULES ---

rule fastqc_se:
   input:
        f1=lambda wildcards: SAMPLES_DICT[wildcards.ID][0],
   output:
        fastqc1=OUTDIR+"/fastqc/{ID}_1_fastqc.zip"
   wildcard_constraints:
        ID="|".join(re.escape(s) for s in SE_SAMPLES) if SE_SAMPLES else "$^"
   shell:
       """
        mkdir -p $(dirname {output.fastqc1})
        
        tmp1=$(mktemp -d)
        fastqc {input.f1} -o $tmp1
        mv $tmp1/*_fastqc.zip {output.fastqc1}
        rm -rf $tmp1
       """

rule fastqcp_se:
   input:
        f1=OUTDIR +"/" +QCTYPE+"/{ID}_1.fq.gz",
   output:
        fastqc1=OUTDIR+"/fastqc/{ID}/1."+QCTYPE+"_fastqc.zip",
   wildcard_constraints:
        ID="|".join(re.escape(s) for s in SE_SAMPLES) if SE_SAMPLES else "$^"
   shell:
       """
        mkdir -p $(dirname {output.fastqc1})
        
        tmp1=$(mktemp -d)
        fastqc {input.f1} -o $tmp1
        mv $tmp1/*_fastqc.zip {output.fastqc1}
        rm -rf $tmp1
       """

rule fastp_se:
   input:
        f1=lambda wildcards: SAMPLES_DICT[wildcards.ID][0],
   output:
        fastqc1=OUTDIR+"/"+QCTYPE+"/{ID}_1.fq",
   params:
        adapters=config.get("adapters", "adapters.fa"),
   log:
        stdout=OUTDIR +"/"+QCTYPE+ "/{ID}/log/stdout.txt",
        stderr=OUTDIR +"/"+QCTYPE+ "/{ID}/log/stderr.txt",
   threads: 8
   wildcard_constraints:
        ID="|".join(re.escape(s) for s in SE_SAMPLES) if SE_SAMPLES else "$^"
   shell:
       """
       mkdir -p $(dirname {output.fastqc1})/log

       if [ "{QCTYPE}" == "fastp" ];
       then
               if [ "{params.adapters}" != "" ] && [ -f "{params.adapters}" ];
               then
                       ADAPTER="--adapter_fasta {params.adapters}" 
               else
                       ADAPTER=""
               fi
               
               fastp --overrepresentation_analysis $ADAPTER --trim_poly_g \
                  --cut_mean_quality 20 --trim_poly_x --length_required 30 --thread {threads} \
                  --cut_window_size 4 --cut_front --cut_tail  \
                  --html $(dirname {output.fastqc1})/{QCTYPE}.html  \
                  --json $(dirname {output.fastqc1})/{QCTYPE}.json \
                  -i {input.f1} -o {output.fastqc1}.tmp > {log.stdout}  2>{log.stderr}

               cutadapt -a 'A{{10}}' -a 'T{{10}}' \
                  -o $(dirname {output.fastqc1})/fastp1c \
                  {output.fastqc1}.tmp >> {log.stdout} 2>>{log.stderr}
                  
               rm {output.fastqc1}.tmp
               
               /mnt/lustre/RDS-live/downing/FASTQC/fastq_quality_trimmer -Q 33 -t 30 -l {MINLEN} -i $(dirname {output.fastqc1})/fastp1c -o {output.fastqc1} -v >> {log.stdout} 2>>{log.stderr}
               
               rm $(dirname {output.fastqc1})/fastp1c
       else
               /mnt/lustre/RDS-live/downing/FASTQC/fastq_quality_trimmer -Q 33 -t 30 -l {MINLEN} -i <(zcat {input.f1}) -o {output.fastqc1} -v > {log.stdout} 2>{log.stderr}
       fi
       """

rule gzip_se:
   input:
        fastqc1=OUTDIR + "/"+QCTYPE+"/{ID}_1.fq",
   output:
        gz1=OUTDIR + "/"+QCTYPE+"/{ID}_1.fq.gz",
   wildcard_constraints:
        ID="|".join(re.escape(s) for s in SE_SAMPLES) if SE_SAMPLES else "$^"
   shell:
       """
       gzip {input.fastqc1}
       """

rule link_se:
   input:
        fastqc1=OUTDIR+"/fastqc/{ID}_1_fastqc.zip",
        f1=OUTDIR+"/fastqc/{ID}/1."+QCTYPE+"_fastqc.zip",
   output:
        OGln1=OUTDIR+"/fastqc/OG/{ID}/1_fastqc.zip",
        FPln1=OUTDIR+"/fastqc/"+QCTYPE+"/{ID}/1_fastqc.zip",
   wildcard_constraints:
        ID="|".join(re.escape(s) for s in SE_SAMPLES) if SE_SAMPLES else "$^"
   shell:
        """
        mkdir -p $(dirname {output.OGln1})
        mkdir -p $(dirname {output.FPln1})
        ln -sf $(readlink -f {input.fastqc1}) {output.OGln1}
        ln -sf $(readlink -f {input.f1}) {output.FPln1}
        """
