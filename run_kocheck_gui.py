#!/usr/bin/env python3
"""
KOCheck GUI Wrapper
A simple graphical interface for running the KOCheck Nextflow pipeline
"""

import tkinter as tk
from tkinter import filedialog, messagebox, ttk
import subprocess
import os
import sys
from pathlib import Path
import threading
import shutil

def check_dependencies():
    """Check for required dependencies (Nextflow, Docker or Conda)"""
    issues = []
    
    # Check for Nextflow
    if shutil.which("nextflow") is None:
        issues.append("Nextflow is not installed. Install from: https://www.nextflow.io/")
    
    # Check for Docker or Conda
    has_docker = shutil.which("docker") is not None
    has_conda = shutil.which("conda") is not None
    
    if not has_docker and not has_conda:
        issues.append("Neither Docker nor Conda is installed.\nPlease install one:")
        issues.append("  • Docker: https://www.docker.com/get-started")
        issues.append("  • Conda: https://docs.conda.io/en/latest/miniconda.html")
    
    return issues

class KOCheckGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("KOCheck - Knockout Validation Pipeline")
        self.root.geometry("700x800")
        
        # Variables to store user inputs
        self.reads_pattern = tk.StringVar()
        self.reference_file = tk.StringVar()
        self.gene_bed_file = tk.StringVar()
        self.marker_fasta_file = tk.StringVar()
        self.marker_contig = tk.StringVar(value="kanR")
        self.output_dir = tk.StringVar(value="results")
        self.delete_ratio = tk.StringVar(value="0.1")
        self.intact_ratio = tk.StringVar(value="0.5")
        self.flank = tk.StringVar(value="500")
        self.min_mapq = tk.StringVar(value="20")
        
        self.setup_ui()
        
    def setup_ui(self):
        """Create the GUI layout"""
        # Title
        title_label = ttk.Label(self.root, text="KOCheck Pipeline", font=("Arial", 16, "bold"))
        title_label.pack(pady=10)
        
        # Create main frame with scrollbar
        main_frame = ttk.Frame(self.root)
        main_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)
        
        # Canvas for scrolling
        canvas = tk.Canvas(main_frame)
        scrollbar = ttk.Scrollbar(main_frame, orient="vertical", command=canvas.yview)
        scrollable_frame = ttk.Frame(canvas)
        
        scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )
        
        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        
        # Required Parameters Section
        required_frame = ttk.LabelFrame(scrollable_frame, text="Required Parameters", padding=10)
        required_frame.pack(fill=tk.X, pady=10)
        
        # Reads pattern
        ttk.Label(required_frame, text="FASTQ Files Pattern:").grid(row=0, column=0, sticky=tk.W, pady=5)
        ttk.Entry(required_frame, textvariable=self.reads_pattern, width=50).grid(row=0, column=1, padx=5)
        ttk.Button(required_frame, text="Browse", command=lambda: self.browse_pattern("reads")).grid(row=0, column=2)
        ttk.Label(required_frame, text="e.g., data/*_{R1,R2}.fq.gz", font=("Arial", 8, "italic")).grid(row=1, column=1, sticky=tk.W)
        
        # Reference genome
        ttk.Label(required_frame, text="Reference Genome:").grid(row=2, column=0, sticky=tk.W, pady=5)
        ttk.Entry(required_frame, textvariable=self.reference_file, width=50).grid(row=2, column=1, padx=5)
        ttk.Button(required_frame, text="Browse", command=lambda: self.browse_file("reference")).grid(row=2, column=2)
        
        # Gene BED file
        ttk.Label(required_frame, text="Target Gene BED:").grid(row=3, column=0, sticky=tk.W, pady=5)
        ttk.Entry(required_frame, textvariable=self.gene_bed_file, width=50).grid(row=3, column=1, padx=5)
        ttk.Button(required_frame, text="Browse", command=lambda: self.browse_file("bed")).grid(row=3, column=2)
        ttk.Label(required_frame, text="(chrom start end format)", font=("Arial", 8, "italic")).grid(row=4, column=1, sticky=tk.W)
        
        # Marker FASTA file
        ttk.Label(required_frame, text="Marker Sequence:").grid(row=5, column=0, sticky=tk.W, pady=5)
        ttk.Entry(required_frame, textvariable=self.marker_fasta_file, width=50).grid(row=5, column=1, padx=5)
        ttk.Button(required_frame, text="Browse", command=lambda: self.browse_file("fasta")).grid(row=5, column=2)
        
        # Optional Parameters Section
        optional_frame = ttk.LabelFrame(scrollable_frame, text="Optional Parameters", padding=10)
        optional_frame.pack(fill=tk.X, pady=10)
        
        # Marker contig name
        ttk.Label(optional_frame, text="Marker Contig Name:").grid(row=0, column=0, sticky=tk.W, pady=5)
        ttk.Entry(optional_frame, textvariable=self.marker_contig, width=30).grid(row=0, column=1, sticky=tk.W, padx=5)
        ttk.Label(optional_frame, text="(default: kanR)", font=("Arial", 8, "italic")).grid(row=0, column=2, sticky=tk.W)
        
        # Output directory
        ttk.Label(optional_frame, text="Output Directory:").grid(row=1, column=0, sticky=tk.W, pady=5)
        ttk.Entry(optional_frame, textvariable=self.output_dir, width=30).grid(row=1, column=1, sticky=tk.W, padx=5)
        ttk.Button(optional_frame, text="Browse", command=lambda: self.browse_directory()).grid(row=1, column=2)
        
        # Thresholds
        ttk.Label(optional_frame, text="Deletion Ratio Threshold:").grid(row=2, column=0, sticky=tk.W, pady=5)
        ttk.Entry(optional_frame, textvariable=self.delete_ratio, width=30).grid(row=2, column=1, sticky=tk.W, padx=5)
        ttk.Label(optional_frame, text="(default: 0.1)", font=("Arial", 8, "italic")).grid(row=2, column=2, sticky=tk.W)
        
        ttk.Label(optional_frame, text="Intact Ratio Threshold:").grid(row=3, column=0, sticky=tk.W, pady=5)
        ttk.Entry(optional_frame, textvariable=self.intact_ratio, width=30).grid(row=3, column=1, sticky=tk.W, padx=5)
        ttk.Label(optional_frame, text="(default: 0.5)", font=("Arial", 8, "italic")).grid(row=3, column=2, sticky=tk.W)
        
        ttk.Label(optional_frame, text="Flank Region (bp):").grid(row=4, column=0, sticky=tk.W, pady=5)
        ttk.Entry(optional_frame, textvariable=self.flank, width=30).grid(row=4, column=1, sticky=tk.W, padx=5)
        ttk.Label(optional_frame, text="(default: 500)", font=("Arial", 8, "italic")).grid(row=4, column=2, sticky=tk.W)
        
        ttk.Label(optional_frame, text="Minimum Mapping Quality:").grid(row=5, column=0, sticky=tk.W, pady=5)
        ttk.Entry(optional_frame, textvariable=self.min_mapq, width=30).grid(row=5, column=1, sticky=tk.W, padx=5)
        ttk.Label(optional_frame, text="(default: 20)", font=("Arial", 8, "italic")).grid(row=5, column=2, sticky=tk.W)
        
        # Pack scrollbar and canvas
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        
        # Buttons
        button_frame = ttk.Frame(self.root)
        button_frame.pack(pady=10, fill=tk.X, padx=10)
        
        ttk.Button(button_frame, text="Run Pipeline", command=self.run_pipeline).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Clear", command=self.clear_fields).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Exit", command=self.root.quit).pack(side=tk.LEFT, padx=5)
        
        # Status/Output frame
        self.output_text = tk.Text(self.root, height=8, width=80)
        self.output_text.pack(pady=5, padx=10, fill=tk.BOTH, expand=True)
        
        # Scrollbar for output text
        scrollbar_text = ttk.Scrollbar(self.output_text)
        scrollbar_text.pack(side=tk.RIGHT, fill=tk.Y)
        self.output_text.config(yscrollcommand=scrollbar_text.set)
        
    def browse_file(self, file_type):
        """Browse for a file"""
        filetypes = {
            "reference": (("FASTA files", "*.fasta *.fa"), ("All files", "*.*")),
            "bed": (("BED files", "*.bed"), ("All files", "*.*")),
            "fasta": (("FASTA files", "*.fasta *.fa"), ("All files", "*.*"))
        }
        
        filename = filedialog.askopenfilename(filetypes=filetypes.get(file_type, []))
        
        if filename:
            if file_type == "reference":
                self.reference_file.set(filename)
            elif file_type == "bed":
                self.gene_bed_file.set(filename)
            elif file_type == "fasta":
                self.marker_fasta_file.set(filename)
    
    def browse_pattern(self, pattern_type):
        """Browse for a directory and create pattern"""
        directory = filedialog.askdirectory()
        if directory:
            pattern = os.path.join(directory, "*_{R1,R2}.fq.gz")
            if pattern_type == "reads":
                self.reads_pattern.set(pattern)
    
    def browse_directory(self):
        """Browse for output directory"""
        directory = filedialog.askdirectory()
        if directory:
            self.output_dir.set(directory)
    
    def clear_fields(self):
        """Clear all input fields"""
        self.reads_pattern.set("")
        self.reference_file.set("")
        self.gene_bed_file.set("")
        self.marker_fasta_file.set("")
        self.marker_contig.set("kanR")
        self.output_dir.set("results")
        self.delete_ratio.set("0.1")
        self.intact_ratio.set("0.5")
        self.flank.set("500")
        self.min_mapq.set("20")
        self.output_text.delete(1.0, tk.END)
    
    def validate_inputs(self):
        """Validate that required fields are filled"""
        if not self.reads_pattern.get():
            messagebox.showerror("Input Error", "Please provide FASTQ files pattern")
            return False
        if not self.reference_file.get():
            messagebox.showerror("Input Error", "Please provide reference genome file")
            return False
        if not self.gene_bed_file.get():
            messagebox.showerror("Input Error", "Please provide target gene BED file")
            return False
        if not self.marker_fasta_file.get():
            messagebox.showerror("Input Error", "Please provide marker sequence file")
            return False
        return True
    
    def run_pipeline(self):
        """Run the Nextflow pipeline"""
        if not self.validate_inputs():
            return
        
        # Check dependencies
        dep_issues = check_dependencies()
        if dep_issues:
            error_msg = "Missing dependencies:\n\n" + "\n".join(dep_issues)
            messagebox.showerror("Missing Dependencies", error_msg)
            return
        
        # Build the nextflow command
        cmd = [
            "nextflow", "run", "main.nf",
            f"--reads '{self.reads_pattern.get()}'",
            f"--reference {self.reference_file.get()}",
            f"--gene_bed {self.gene_bed_file.get()}",
            f"--marker_fasta {self.marker_fasta_file.get()}",
            f"--marker_contig {self.marker_contig.get()}",
            f"--outdir {self.output_dir.get()}",
            f"--delete_ratio {self.delete_ratio.get()}",
            f"--intact_ratio {self.intact_ratio.get()}",
            f"--flank {self.flank.get()}",
            f"--min_mapq {self.min_mapq.get()}"
        ]
        
        # Run in a separate thread to prevent GUI freezing
        thread = threading.Thread(target=self.execute_pipeline, args=(cmd,))
        thread.daemon = True
        thread.start()
    
    def execute_pipeline(self, cmd):
        """Execute the pipeline command"""
        self.output_text.insert(tk.END, "Starting KOCheck pipeline...\n")
        self.output_text.insert(tk.END, f"Command: {' '.join(cmd)}\n\n")
        self.output_text.see(tk.END)
        self.root.update()
        
        try:
            # Join command into a single string for shell execution
            cmd_str = " ".join(cmd)
            process = subprocess.Popen(
                cmd_str,
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                universal_newlines=True,
                bufsize=1
            )
            
            # Stream output in real-time
            for line in iter(process.stdout.readline, ""):
                if line:
                    self.output_text.insert(tk.END, line)
                    self.output_text.see(tk.END)
                    self.root.update()
            
            process.wait()
            
            if process.returncode == 0:
                self.output_text.insert(tk.END, "\n✓ Pipeline completed successfully!\n")
                messagebox.showinfo("Success", "Pipeline completed successfully!")
            else:
                self.output_text.insert(tk.END, f"\n✗ Pipeline failed with error code {process.returncode}\n")
                messagebox.showerror("Error", f"Pipeline failed with error code {process.returncode}")
                
        except Exception as e:
            error_msg = f"Error running pipeline: {str(e)}\n"
            self.output_text.insert(tk.END, error_msg)
            messagebox.showerror("Error", error_msg)
        
        self.output_text.see(tk.END)

if __name__ == "__main__":
    root = tk.Tk()
    
    # Check for dependencies on startup
    dep_issues = check_dependencies()
    if dep_issues:
        root.withdraw()  # Hide the main window initially
        error_msg = "Missing required dependencies:\n\n" + "\n".join(dep_issues) + "\n\nPlease install the missing tools and try again."
        messagebox.showerror("Setup Required", error_msg)
        root.destroy()
        sys.exit(1)
    
    app = KOCheckGUI(root)
    root.mainloop()
