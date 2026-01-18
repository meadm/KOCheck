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

class ToolTip:
    """Create a tooltip for a given widget"""
    def __init__(self, widget, text):
        self.widget = widget
        self.text = text
        self.tooltip = None
        widget.bind("<Enter>", self.enter)
        widget.bind("<Leave>", self.leave)
    
    def enter(self, event=None):
        if self.tooltip:
            return
        x = self.widget.winfo_rootx() + self.widget.winfo_width()
        y = self.widget.winfo_rooty()
        self.tooltip = tk.Toplevel(self.widget)
        self.tooltip.wm_overrideredirect(True)
        self.tooltip.wm_geometry(f"+{x}+{y}")
        label = tk.Label(self.tooltip, text=self.text, background="#ffffe0", 
                        relief=tk.SOLID, borderwidth=1, font=("Arial", 10))
        label.pack()
    
    def leave(self, event=None):
        if self.tooltip:
            self.tooltip.destroy()
            self.tooltip = None

class PlaceholderEntry(ttk.Entry):
    """Entry widget with placeholder text that disappears on focus"""
    def __init__(self, parent, placeholder="", *args, **kwargs):
        super().__init__(parent, *args, **kwargs)
        self.placeholder = placeholder
        self.placeholder_active = False
        self.default_color = self["foreground"]
        
        self.bind("<FocusIn>", self._on_focus_in)
        self.bind("<FocusOut>", self._on_focus_out)
        
        self._show_placeholder()
    
    def _show_placeholder(self):
        """Display placeholder text"""
        if not self.get():
            self.placeholder_active = True
            self.insert(0, self.placeholder)
            self.config(foreground="gray")
    
    def _on_focus_in(self, event=None):
        """Remove placeholder when entry gains focus"""
        if self.placeholder_active:
            self.delete(0, tk.END)
            self.config(foreground=self.default_color)
            self.placeholder_active = False
    
    def _on_focus_out(self, event=None):
        """Restore placeholder if entry is empty"""
        if not self.get():
            self._show_placeholder()

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
        self.resume = tk.BooleanVar(value=False)
        
        self.setup_ui()
        
    def setup_ui(self):
        """Create the GUI layout"""
        # Title
        title_label = ttk.Label(self.root, text="KOCheck ✅", font=("Arial", 16, "bold"))
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
        reads_label = ttk.Label(required_frame, text="FASTQ Files Pattern:")
        reads_label.grid(row=0, column=0, sticky=tk.W, pady=5)
        ToolTip(reads_label, "Glob pattern to your paired-end FASTQ files")
        PlaceholderEntry(required_frame, placeholder="path/to/*_{R1,R2}.fq.gz", textvariable=self.reads_pattern, width=50).grid(row=0, column=1, padx=5)
        ttk.Button(required_frame, text="Browse", command=lambda: self.browse_pattern("reads")).grid(row=0, column=2)
        
        # Reference genome
        ref_label = ttk.Label(required_frame, text="Reference Genome:")
        ref_label.grid(row=2, column=0, sticky=tk.W, pady=5)
        ToolTip(ref_label, "FASTA file of reference genome used as alignment target")
        PlaceholderEntry(required_frame, placeholder="path/to/reference.fasta", textvariable=self.reference_file, width=50).grid(row=2, column=1, padx=5)
        ttk.Button(required_frame, text="Browse", command=lambda: self.browse_file("reference")).grid(row=2, column=2)
        
        # Gene BED file
        bed_label = ttk.Label(required_frame, text="Target Gene BED:")
        bed_label.grid(row=3, column=0, sticky=tk.W, pady=5)
        ToolTip(bed_label, "Single-line BED file with target gene coordinates (format: chromosome start end)")
        PlaceholderEntry(required_frame, placeholder="path/to/target_gene.bed", textvariable=self.gene_bed_file, width=50).grid(row=3, column=1, padx=5)
        ttk.Button(required_frame, text="Browse", command=lambda: self.browse_file("bed")).grid(row=3, column=2)
        #ttk.Label(required_frame, text="(format: chromosome start end)", font=("Arial", 10, "italic")).grid(row=4, column=1, sticky=tk.W)
        
        # Marker FASTA file
        marker_label = ttk.Label(required_frame, text="Marker Sequence:")
        marker_label.grid(row=5, column=0, sticky=tk.W, pady=5)
        ToolTip(marker_label, "FASTA file with resistance marker sequence")
        PlaceholderEntry(required_frame, placeholder="path/to/marker.fasta", textvariable=self.marker_fasta_file, width=50).grid(row=5, column=1, padx=5)
        ttk.Button(required_frame, text="Browse", command=lambda: self.browse_file("fasta")).grid(row=5, column=2)
        
        # Optional Parameters Section
        optional_frame = ttk.LabelFrame(scrollable_frame, text="Optional Parameters", padding=10)
        optional_frame.pack(fill=tk.X, pady=10)
        
        # Resume previous run (moved to top)
        resume_check = ttk.Checkbutton(optional_frame, text="Resume from previous run (use cached results)", 
                        variable=self.resume)
        resume_check.grid(row=0, column=0, columnspan=3, sticky=tk.W, pady=10)
        ToolTip(resume_check, "Enable to use Nextflow's caching and skip previously completed steps")
        
        # Marker contig name
        marker_name_label = ttk.Label(optional_frame, text="Marker Contig Name:")
        marker_name_label.grid(row=1, column=0, sticky=tk.W, pady=5)
        ToolTip(marker_name_label, "Name/ID of marker contig in the marker file (default: kanR)")
        PlaceholderEntry(optional_frame, placeholder="kanR", textvariable=self.marker_contig, width=30).grid(row=1, column=1, sticky=tk.W, padx=5)
        ttk.Label(optional_frame, text="(default: kanR)", font=("Arial", 12, "italic")).grid(row=1, column=2, sticky=tk.W)
        
        # Output directory
        output_label = ttk.Label(optional_frame, text="Output Directory:")
        output_label.grid(row=2, column=0, sticky=tk.W, pady=5)
        ToolTip(output_label, "Directory where pipeline results will be saved")
        PlaceholderEntry(optional_frame, placeholder="./results", textvariable=self.output_dir, width=30).grid(row=2, column=1, sticky=tk.W, padx=5)
        ttk.Button(optional_frame, text="Browse", command=lambda: self.browse_directory()).grid(row=2, column=2)
        
        # Thresholds
        delete_ratio_label = ttk.Label(optional_frame, text="Deletion Ratio Threshold:")
        delete_ratio_label.grid(row=3, column=0, sticky=tk.W, pady=5)
        ToolTip(delete_ratio_label, "Coverage threshold for deletion classification (default: 0.1)")
        PlaceholderEntry(optional_frame, placeholder="0.1", textvariable=self.delete_ratio, width=30).grid(row=3, column=1, sticky=tk.W, padx=5)
        ttk.Label(optional_frame, text="(default: 0.1)", font=("Arial", 12, "italic")).grid(row=3, column=2, sticky=tk.W)
        
        ttk.Label(optional_frame, text="Intact Ratio Threshold:").grid(row=4, column=0, sticky=tk.W, pady=5)
        PlaceholderEntry(optional_frame, placeholder="0.5", textvariable=self.intact_ratio, width=30).grid(row=4, column=1, sticky=tk.W, padx=5)
        ttk.Label(optional_frame, text="(default: 0.5)", font=("Arial", 12, "italic")).grid(row=4, column=2, sticky=tk.W)
        intact_ratio_label = ttk.Label(optional_frame, text="Intact Ratio Threshold:")
        ToolTip(intact_ratio_label, "Coverage threshold for intact gene classification (default: 0.5)")
        
        ttk.Label(optional_frame, text="Flank Region (bp):").grid(row=5, column=0, sticky=tk.W, pady=5)
        PlaceholderEntry(optional_frame, placeholder="500", textvariable=self.flank, width=30).grid(row=5, column=1, sticky=tk.W, padx=5)
        ttk.Label(optional_frame, text="(default: 500)", font=("Arial", 12, "italic")).grid(row=5, column=2, sticky=tk.W)
        flank_label = ttk.Label(optional_frame, text="Flank Region (bp):")
        ToolTip(flank_label, "Flanking region in bp for visualization and analysis (default: 500)")
        
        ttk.Label(optional_frame, text="Minimum Mapping Quality:").grid(row=6, column=0, sticky=tk.W, pady=5)
        PlaceholderEntry(optional_frame, placeholder="20", textvariable=self.min_mapq, width=30).grid(row=6, column=1, sticky=tk.W, padx=5)
        ttk.Label(optional_frame, text="(default: 20)", font=("Arial", 12, "italic")).grid(row=6, column=2, sticky=tk.W)
        mapq_label = ttk.Label(optional_frame, text="Minimum Mapping Quality:")
        ToolTip(mapq_label, "Minimum mapping quality for junction analysis (default: 20)")
        
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
        self.resume.set(False)
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
        
        # Add -resume flag if checked
        if self.resume.get():
            cmd.insert(2, "-resume")
        
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
                self.output_text.insert(tk.END, "\n" + "="*60 + "\n")
                self.output_text.insert(tk.END, "✓ Pipeline completed successfully!\n")
                self.output_text.insert(tk.END, "="*60 + "\n")
                self.output_text.see(tk.END)
                self.root.update()
                messagebox.showinfo("Success", "Pipeline completed successfully!\n\nResults are ready in the output directory.")
            else:
                self.output_text.insert(tk.END, "\n" + "="*60 + "\n")
                self.output_text.insert(tk.END, f"✗ Pipeline failed with error code {process.returncode}\n")
                self.output_text.insert(tk.END, "="*60 + "\n")
                self.output_text.see(tk.END)
                self.root.update()
                messagebox.showerror("Pipeline Error", f"Pipeline failed with error code {process.returncode}\n\nCheck the output above for details.")
                
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
