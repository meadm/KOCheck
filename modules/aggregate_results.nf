process AGGREGATE_RESULTS {
    tag "aggregate_results"
    
    publishDir "${params.outdir}/summary", mode: 'copy'

    input:
    path all_files

    output:
    path "kocheck_summary.csv", emit: summary_csv
    path "kocheck_report.html", emit: report_html

    script:
    """
    set -euo pipefail

    python3 << 'PYTHON_EOF'
import pandas as pd
import base64
import os
import sys
import glob
from pathlib import Path

# Read all marker_check module output CSV files
csv_files = sorted(glob.glob("*.marker_check.csv"))
if not csv_files:
    print("ERROR: No marker_check CSV files found", file=sys.stderr)
    sys.exit(1)

# Combine all CSV files
dfs = []
for csv_file in csv_files:
    df = pd.read_csv(csv_file)
    dfs.append(df)

combined_df = pd.concat(dfs, ignore_index=True)

# Sort by sample_id
combined_df = combined_df.sort_values('sample_id').reset_index(drop=True)

# Reorder columns: sample_id, status, deletion_status, marker_present, marker_ratio, marker_copy, junction_support
column_order = ['sample_id', 'status', 'deletion_status', 'marker_present', 'marker_ratio', 'marker_copy', 'junction_support']
combined_df = combined_df[column_order]

# Write combined CSV
combined_df.to_csv("kocheck_summary.csv", index=False)
print(f"Created summary CSV with {len(combined_df)} samples")

# Create HTML report
html_content = []
html_content.append('''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>KOCheck Summary Report</title>
    <style>
        body {
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, sans-serif;
            margin: 0;
            padding: 20px;
            background-color: #f5f5f5;
        }
        .container {
            max-width: 1400px;
            margin: 0 auto;
            background-color: white;
            padding: 30px;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        h1 {
            color: #2c3e50;
            border-bottom: 3px solid #3498db;
            padding-bottom: 10px;
        }
        table {
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
            font-size: 14px;
        }
        th {
            background-color: #3498db;
            color: white;
            padding: 12px;
            text-align: left;
            font-weight: 600;
        }
        td {
            padding: 10px 12px;
            border-bottom: 1px solid #ddd;
        }
        tr:hover {
            background-color: #f8f9fa;
        }
        .status-PASS {
            background-color: #d4edda;
            color: #155724;
            font-weight: bold;
        }
        .status-FAIL {
            background-color: #f8d7da;
            color: #721c24;
            font-weight: bold;
        }
        .status-REVIEW {
            background-color: #fff3cd;
            color: #856404;
            font-weight: bold;
        }
        .deleted {
            color: #28a745;
            font-weight: bold;
        }
        .intact {
            color: #dc3545;
            font-weight: bold;
        }
        .ambiguous {
            color: #6c757d;
        }
        .plot-container {
            margin: 30px 0;
            padding: 20px;
            background-color: #f8f9fa;
            border-radius: 4px;
        }
        .plot-container h2 {
            color: #2c3e50;
            margin-top: 0;
        }
        .plot-container img {
            max-width: 100%;
            height: auto;
            border: 1px solid #ddd;
            border-radius: 4px;
        }
        .summary-stats {
            display: flex;
            gap: 20px;
            margin: 20px 0;
            flex-wrap: wrap;
        }
        .stat-box {
            flex: 1;
            min-width: 150px;
            padding: 15px;
            background-color: #f8f9fa;
            border-radius: 4px;
            border-left: 4px solid #3498db;
        }
        .stat-box h3 {
            margin: 0 0 5px 0;
            color: #2c3e50;
            font-size: 14px;
        }
        .stat-box p {
            margin: 0;
            font-size: 24px;
            font-weight: bold;
            color: #3498db;
        }
    </style>
</head>
<body>
    <div class="container">
        <h1>KOCheck Summary Report</h1>
''')

# Add summary statistics
total_samples = len(combined_df)
status_counts = combined_df['status'].value_counts()
deletion_counts = combined_df['deletion_status'].value_counts()

# Define the order for status display
status_order = ['PASS', 'REVIEW', 'FAIL']

html_content.append('''        <div class="summary-stats">
            <div class="stat-box">
                <h3>Total Samples</h3>
                <p>{}</p>
            </div>
'''.format(total_samples))

# Display status counts in specified order
for status in status_order:
    if status in status_counts:
        count = status_counts[status]
        html_content.append('''            <div class="stat-box">
                <h3>{}</h3>
                <p>{}</p>
            </div>
'''.format(status, count))

html_content.append('''        </div>
''')

# Add summary table
html_content.append('''        <h2>Sample Summary</h2>
        <table>
            <thead>
                <tr>
                    <th>Sample ID</th>
                    <th>Overall Status</th>
                    <th>Deletion Status</th>
                    <th>Marker Present</th>
                    <th>Marker Ratio</th>
                    <th>Marker Copy</th>
                    <th>Junction Support</th>
                </tr>
            </thead>
            <tbody>
''')

for _, row in combined_df.iterrows():
    sample_id = row['sample_id']
    deletion_status = row['deletion_status']
    marker_present = row['marker_present']
    marker_copy = row['marker_copy']
    marker_ratio = f"{float(row['marker_ratio']):.3f}"
    junction_support = int(row['junction_support'])
    status = row['status']
    
    # Determine CSS classes
    deletion_class = deletion_status.lower()
    status_class = f"status-{status}"
    
    html_content.append(f'''                <tr>
                    <td><strong>{sample_id}</strong></td>
                    <td class="{status_class}">{status}</td>
                    <td class="{deletion_class}">{deletion_status}</td>
                    <td>{marker_present}</td>
                    <td>{marker_ratio}</td>
                    <td>{marker_copy}</td>
                    <td>{junction_support}</td>
                </tr>
''')

html_content.append('''            </tbody>
        </table>
''')

# Add coverage plots
plot_files = sorted(glob.glob("*.gene_coverage.png"))
if plot_files:
    html_content.append('''        <h2>Coverage Plots</h2>
''')
    
    for plot_file in plot_files:
        sample_id = plot_file.replace('.gene_coverage.png', '')
        
        # Read and encode image
        with open(plot_file, 'rb') as f:
            img_data = base64.b64encode(f.read()).decode('utf-8')
        
        html_content.append(f'''        <div class="plot-container">
            <h2>{sample_id}</h2>
            <img src="data:image/png;base64,{img_data}" alt="Coverage plot for {sample_id}">
        </div>
''')

html_content.append('''    </div>
</body>
</html>
''')

# Write HTML file
with open("kocheck_report.html", "w") as f:
    f.write("".join(html_content))

print(f"Created HTML report with {len(plot_files)} coverage plots")
PYTHON_EOF
    """
}

