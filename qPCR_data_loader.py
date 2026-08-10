import os
import re
import pandas as pd
import importlib

try:
    pd = importlib.import_module("pandas")
except ImportError as exc:
    raise ImportError("pandas is required to run qPCR_data_loader.py") from exc


def load_sample_name_mapping(raw_filepath):
    """
    Load the first sample name seen for each plate position from a raw qPCR file.
    """
    position_map = {}
    with open(raw_filepath, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('Raw Data') or line.startswith('SamplePos'):
                continue
            parts = line.split()
            if len(parts) < 8:
                continue
            sample_pos = parts[0]
            if sample_pos in position_map:
                continue
            sample_name = " ".join(parts[1:-6])
            if sample_name == "":
                sample_name = " ".join(parts[1:-5])
            position_map[sample_pos] = sample_name
    return position_map


def load_qpcr_data(filepath, raw_sample_file=None):
    """
    Load qPCR data from the summarized results file.
    Handles irregular whitespace formatting and optional sample name mapping.
    
    Parameters:
    -----------
    filepath : str
        Path to the qPCR results file (e.g., "230913 qPCR and SD.txt")
    raw_sample_file : str, optional
        Path to the raw qPCR file containing sample names by well position
    
    Returns:
    --------
    df : pandas.DataFrame
        DataFrame with qPCR data including sample information and Ct values
    """
    data = []
    with open(filepath, 'r', encoding='utf-8') as f:
        lines = f.readlines()
        for line in lines[1:]:
            fields = [field.strip() for field in line.split('\t') if field.strip()]
            if len(fields) >= 3:
                try:
                    samples = fields[0]
                    mean_cp = float(fields[1])
                    std_cp = float(fields[2])
                    mean_conc = float(fields[3]) if len(fields) > 3 else 0
                    std_conc = float(fields[4]) if len(fields) > 4 else 0
                    data.append({
                        'Samples': samples,
                        'MeanCp': mean_cp,
                        'STD_Cp': std_cp,
                        'Mean_conc': mean_conc,
                        'STD_conc': std_conc
                    })
                except ValueError:
                    continue
    df = pd.DataFrame(data)
    if len(df) > 0:
        df[['Well1', 'Well2']] = df['Samples'].str.split(', ', expand=True)
        df['Well1'] = df['Well1'].str.strip()
        df['Well2'] = df['Well2'].str.strip()
        df['Sample_Number'] = df['Well1'].str.extract(r'(\d+)$')[0].astype(int)
        if raw_sample_file:
            mapping = load_sample_name_mapping(raw_sample_file)
            df['Sample_Name'] = df['Well1'].map(mapping)
            df['Sample_Name'] = df['Sample_Name'].fillna(df['Well2'].map(mapping))
            df['Sample_Name'] = df['Sample_Name'].fillna(df['Samples'])
        else:
            df['Sample_Name'] = df['Samples']
        df['Day'] = df['Sample_Name'].str.extract(r'^(D[12])', expand=False).fillna('Unknown')
    return df


def parse_raw_plate_data(filepath):
    """Parse the raw export from the qPCR instrument into well-level metadata."""
    if not filepath or not os.path.exists(filepath):
        return []

    records = []
    with open(filepath, 'r', encoding='utf-8') as f:
        for line in f:
            if not line.strip() or line.startswith('Experiment:') or line.startswith('Include'):
                continue
            fields = [field.strip() for field in line.split('\t')]
            if len(fields) < 7:
                continue

            include = fields[0].lower() == 'true'
            if not include:
                continue

            position = fields[2]
            name = fields[3]
            cp = fields[4]
            concentration = fields[5] if fields[5] else None
            status = fields[6] if len(fields) > 6 else ''
            day_match = re.search(r'(D\d+)', name, flags=re.IGNORECASE)
            day = day_match.group(1).upper() if day_match else 'Unknown'
            day_number = int(re.search(r'(\d+)', day).group(1)) if day.startswith('D') else None

            try:
                cp_value = float(cp) if cp else None
            except ValueError:
                cp_value = None

            try:
                conc_value = float(concentration) if concentration else None
            except ValueError:
                conc_value = None

            records.append({
                'position': position,
                'well': position,
                'name': name,
                'sample_name': name,
                'cp': cp_value,
                'concentration': conc_value,
                'day': day,
                'day_number': day_number,
                'status': status,
            })

    return records


def separate_genes(df):
    """
    Separate target and reference genes based on well positions.
    
    Assumes:
    - Target gene: Wells A and B rows
    - Reference gene: Wells C and D rows
    
    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame from load_qpcr_data()
    
    Returns:
    --------
    target_df : pandas.DataFrame
        Data for target gene (rows A/B)
    reference_df : pandas.DataFrame
        Data for reference gene (rows C/D)
    """
    # Extract row letter from well position
    df['Row'] = df['Well1'].str[0]
    
    # Separate by row
    target_df = df[df['Row'].isin(['A', 'B'])].copy()
    reference_df = df[df['Row'].isin(['C', 'D'])].copy()
    
    # Reset indices
    target_df = target_df.reset_index(drop=True)
    reference_df = reference_df.reset_index(drop=True)
    
    return target_df, reference_df


def build_sample_key(row):
    """Create a stable identifier for a sample across day and sample number."""
    sample_number = row.get('Sample_Number')
    day = str(row.get('Day') or '').strip()
    sample_name = str(row.get('Sample_Name') or '').strip()

    if pd.notna(sample_number):
        return f"{day if day and day.lower() != 'unknown' else 'unknown'}|{sample_number}"
    if sample_name:
        return sample_name

    well = row.get('Well1') or row.get('Target_Well') or row.get('Reference_Well') or row.get('Sample_Name')
    return str(well or row.get('Samples') or '').strip()


def make_well_raw_lookup(df, label='sample'):
    """Build a per-well raw-data lookup from a loaded qPCR table."""
    if df is None or df.empty:
        return {}

    lookup = {}
    for _, row in df.iterrows():
        well = row.get('Well1') or row.get('Samples')
        if pd.isna(well) or not str(well).strip():
            continue

        sample_key = build_sample_key(row)
        lookup[str(well).strip()] = {
            'label': label,
            'sample_key': sample_key,
            'sample_number': int(row['Sample_Number']) if 'Sample_Number' in row and pd.notna(row.get('Sample_Number')) else None,
            'well': str(well).strip(),
            'sample_name': str(row.get('Sample_Name') or row.get('Samples') or str(well)).strip(),
            'day': str(row.get('Day') or '').strip() or None,
            'mean_cp': float(row.get('MeanCp')) if pd.notna(row.get('MeanCp')) else None,
            'std_cp': float(row.get('STD_Cp')) if pd.notna(row.get('STD_Cp')) else None,
            'mean_conc': float(row.get('Mean_conc')) if pd.notna(row.get('Mean_conc')) else None,
            'std_conc': float(row.get('STD_conc')) if pd.notna(row.get('STD_conc')) else None,
        }

    return lookup


def get_sample_groups(df):
    """
    Extract sample information from well names.
    Creates a mapping of samples to their Ct values.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame from load_qpcr_data()
    
    Returns:
    --------
    groups : dict
        Dictionary with sample information
    """
    # Extract sample group from well position (column number)
    df['Sample_Number'] = df['Well1'].str[1:].astype(int)
    
    groups = {
        'sample_numbers': sorted(df['Sample_Number'].unique()),
        'data': df
    }
    
    return groups


def display_data_summary(df, target_df, reference_df):
    """
    Display a summary of the loaded data.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        Full dataset
    target_df : pandas.DataFrame
        Target gene data
    reference_df : pandas.DataFrame
        Reference gene data
    """
    print("=" * 60)
    print("qPCR DATA SUMMARY")
    print("=" * 60)
    print(f"\nTotal samples: {len(df)}")
    print(f"Target gene samples (A/B): {len(target_df)}")
    print(f"Reference gene samples (C/D): {len(reference_df)}")
    
    print("\n--- TARGET GENE (A/B) ---")
    print(f"Mean Ct range: {target_df['MeanCp'].min():.2f} - {target_df['MeanCp'].max():.2f}")
    print(f"Mean STD Cp: {target_df['STD_Cp'].mean():.4f}")
    
    print("\n--- REFERENCE GENE (C/D) ---")
    print(f"Mean Ct range: {reference_df['MeanCp'].min():.2f} - {reference_df['MeanCp'].max():.2f}")
    print(f"Mean STD Cp: {reference_df['STD_Cp'].mean():.4f}")
    print("\n" + "=" * 60)


if __name__ == "__main__":
    # Example usage
    filepath = "230913 qPCR and SD.txt"
    
    # Load data
    df = load_qpcr_data(filepath)
    print("Raw data loaded:")
    print(df)
    
    # Separate genes
    target_df, reference_df = separate_genes(df)
    
    # Display summary
    display_data_summary(df, target_df, reference_df)
    
    # Get sample groups
    groups = get_sample_groups(df)
    print(f"\nSample numbers: {groups['sample_numbers']}")
